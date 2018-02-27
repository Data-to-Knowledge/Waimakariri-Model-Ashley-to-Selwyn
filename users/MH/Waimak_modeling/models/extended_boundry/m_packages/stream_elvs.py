import numpy as np
import pandas as pd
import geopandas as gpd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
from copy import deepcopy

def clean_weird_points(data):
    """
    take the reach data and clean up any reaches that are uphill and downstream
    :param data: sfr reach data
    :return:
    """
    counter = 0
    while any(np.array(data.iloc[0:-2])<=np.array(data.iloc[1:-1])):
        counter +=1
        if counter > 10:
            raise ValueError('run away while loop')
        for i in data.index:
            if i==data.index.min():
                continue
            elif i==data.index.max():
                continue

            prev = data.loc[i-1]
            if data.loc[i] >= prev:
                next = data.loc[i+1]
                data.loc[i] = (prev+next)/2
    return data

def get_reach_elv():
    """
    gets the stream tops to use for the reaches either from the previous model or lidar
    :return:
    """

    str_seg_trib_dict = { # dictionary to link stream segments to tributaries
        1: None,
        2: (1,),
        3: None,
        4: (2,),
        5: None,
        6: (3,),
        7: None,
        8: None,
        9: (6,),
        10: (7,),
        11: (4,),
        12: (5,8),
        13: (10,),
        14: (9,),
        15: None,
        16: (13,),
        17: (11,),
        18: (12,15),
        19: None,
        20: (14,),
        21: (16,19),
        22: None,
        23: (21,),
        24: None,
        25: None,
        26: (25,43),
        27: None,
        28: None,
        29: (17,20),
        30: (18,22),
        31: (24,30),
        32: (28,31),
        33: (32,),
        34: (33,),
        35: (29,),
        36: (27,44),
        37: (36,),
        38: None,
        39: None,
        40: None,
        41: None,
        42: (23,38,41),
        43: (39,42),
        44: (26,40)}
    down_str_seg_dict = {}
    for key in str_seg_trib_dict.keys():
        if str_seg_trib_dict[key] is None:
            continue

        tribs = str_seg_trib_dict[key]
        for i in tribs:
            down_str_seg_dict[i] = key




    data = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\stream_drn_values\Streams_MH.shp")
    data.loc[:,'NCS_elev'] = data.loc[:,'NCS_elev'].astype(float)
    data.loc[(data.iseg==3)&(data.ireach==94),'NCS_elev'] = 183


    #for new/old I could set

    elv = smt.calc_elv_db()
    for i in data.index:
        row, col = data.loc[i,['i','j']]
        data.loc[i,'ground'] = elv[0,row,col]
    xs, ys = smt.get_model_x_y()
    data.loc[data['NCS_elev']>data.ground,'NCS_elev'] = data.loc[data['NCS_elev']>data.ground,'ground']
    data.loc[data['strtop']>data.ground,'strtop'] = data.loc[data['strtop']>data.ground,'ground']

    data.loc[(data.iseg == 16)& (data.ireach == 1),'strtop'] = 107


    segments = smt.df_to_array(data,'iseg')
    no_flow = smt.get_no_flow(0)
    old='old'
    new = 'new'
    seg_dict = {1:old,
     2:old,
     3:new,
     4:old,
     5:old,
     6:new,
     7:old,
     8:old,
     9:old,
     10:old,
     11:old,
     12:old,
     13:new,
     14:old,
     15:old,
     16:old,
     17:old,
     18:old,
     19:old,
     20:old,
     21:old,
     22:old,
     23:old,
     24:old,
     25:old,
     26:old,
     27:old,
     28:old,
     29:old,
     30:old,
     31:old,
     32:old,
     33:old,
     34:old,
     35:old,
     36:new,
     37:new,
     38:old,
     39:old,
     40:new,
     41:old,
     42:old,
     43:old,
     44:new,
     }

    plot = False # hold over to check stuff
    if plot:
        for seg in set(data.iseg):
            if seg not in [3]:
                continue
            temp = data.loc[data.iseg==seg]
            if seg_dict[seg] =='new':
                temp.loc[:,'elev_to_use'] = clean_weird_points(temp.loc[:,'NCS_elev'])
            elif seg_dict[seg]=='old':
                temp.loc[:,'elev_to_use'] = clean_weird_points(temp.loc[:,'strtop'])
            else:
                raise ValueError('shouldnt get here')


            fig,(ax, ax2)=plt.subplots(1, 2, figsize=(18.5, 9.5))
            ax.plot(temp.ireach,temp.ground,color='k',label='ground')
            ax.plot(temp.ireach,temp.strtop,color='r',label='old_top')
            ax.plot(temp.ireach,temp['NCS_elev'],color='g',label='new_top')
            ax.set_title('segment {}'.format(seg))
            ax.legend()
            ax2.pcolormesh(xs,ys,np.isclose(segments,seg))
            ax2.set_aspect('equal')
            ax2.contour(xs[0:190,:], ys[0:190,:], no_flow[0:190,:])
            plt.show(fig)


    #add seg_to_use to data
    for seg in set(data.iseg):
        temp = data.loc[data.iseg==seg]
        if seg_dict[seg] == 'new':
            temp.loc[:, 'elev_to_use'] = clean_weird_points(temp.loc[:, 'NCS_elev'])
        elif seg_dict[seg] == 'old':
            temp.loc[:, 'elev_to_use'] = clean_weird_points(temp.loc[:, 'strtop'])
        data.loc[data.iseg == seg,'elev_to_use'] = temp.loc[:,'elev_to_use']

    # check for problems across segments
    for seg in set(data.iseg):

        if seg_dict[seg] == 'new':
            temp = data.loc[data.iseg == seg]
            bot = temp.loc[temp.ireach==temp.ireach.max(),'elev_to_use'].iloc[0]
            # check downstream
            try:
                down_seg = down_str_seg_dict[seg]
                down = data.loc[data.iseg == down_seg]
                down_top = down.loc[down.ireach==down.ireach.min(),'elev_to_use'].iloc[0]
                if down_top > bot:
                    raise ValueError('segment: {}, bottom {}, down_seg {} top of downstream {}'.format(seg, bot,
                                                                                                       down_seg,
                                                                                                       down_top))
            except KeyError as val:
                pass


            top = temp.loc[temp.ireach==temp.ireach.min(),'elev_to_use'].iloc[0]
            # check upstream
            up_segs = str_seg_trib_dict[seg]
            if up_segs is None:
                continue
            for us in up_segs:
                up = data.loc[data.iseg == us]
                up_bot = up.loc[up.ireach==up.ireach.max(),'elev_to_use'].iloc[0]
                if up_bot < top:
                    raise ValueError('segment{} top: {}, upseg: {} upbot: {}'.format(seg,top,us,up_bot))

    return data.loc[:,'elev_to_use']


if __name__ == '__main__':
    #tests
    get_reach_elv()