"""
Author: matth
Date Created: 26/04/2017 4:04 PM
"""

from __future__ import division
import pandas as pd
import numpy as np
import flopy_mh as flopy
# where segments overlap cell by cell flow is likely doubled look at further.
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

#carefully look at behavior when streams don't flood

def calc_full_stream_conc(mt3d, stress_period, time_step):
    # return a dataframe will full stream concentrations for all model cells with streams
    # return a pd.Dataframe(index = (layer,row,col), columns = ['segment','reach','flow_gw','flow_sw', 'c_gw', 'c_sw']

    segments = range(1,45)
    out_data = _fill_stream_conc_data(segments, mt3d, stress_period, time_step)
    # assume that there are no subterrainan streams
    #  watch trib overlap
    out_data.index = pd.MultiIndex.from_tuples([(0,out_data.loc[m,'i'],out_data.loc[m,'j']) for m in out_data.index])

    return out_data


def calc_stream_conc_at_pt(row, col, mt3d, stress_period, time_step):
    # as above but only for a single point (not sure which is quicker)
    m = mt3d.mf
    str_data = m.str.stress_period_data.data[stress_period]
    seg = str_data['segment'][(str_data['i']==row) & (str_data['j']==col)]
    if len(seg) > 1:
        raise ValueError('more than one segment identified at row: {}, col: {}. likely at a junction'.format(row,col))
    elif len(seg) == 0:
        raise ValueError('no stream boundry condition at row: {}, col: {}'.format(row,col))

    segments =get_all_segs(seg[0])

    out_data = _fill_stream_conc_data(segments,mt3d,stress_period,time_step)
    # assume that there are no subterrainan streams
    # watch trib overlap
    out_data.index = pd.MultiIndex.from_tuples([(0,out_data.loc[m,'i'],out_data.loc[m,'j']) for m in out_data.index])

    return out_data


def _fill_stream_conc_data (segments, mt3d, stress_period, time_step):
    #function to fill the out array based on this record
    m = mt3d.mf
    # initalize data frame
    str_data = m.str.stress_period_data.data[stress_period]
    out_data = pd.DataFrame(str_data)
    out_data = out_data.drop([ u'cond', u'sbot', u'stop', u'width', u'slope', u'rough'],axis=1)
    out_data['flow_gw','flow_sw', 'c_gw', 'c_sw'] = np.nan
    out_data.index = pd.MultiIndex.from_tuples([(out_data.segment.loc[i],out_data.reach.loc[i])
                                                for i in out_data.index])

    # load data from model
    str_flow_data = flopy.utils.CellBudgetFile('{}/{}'.format(m.model_ws,m.get_output(unit=744)
                                                              )).get_data(kstpkper=(time_step,stress_period),
                                                                          full3D=True)[0]
    str_recharge_data = flopy.utils.CellBudgetFile('{}/{}'.format(m.model_ws,m.get_output(unit=743)
                                                              )).get_data(kstpkper=(time_step,stress_period),
                                                                          full3D=True)[0]
    gw_conc_data = None # once MT3D is up and running

    # load flow and concentration data into the dataframe
    for i in out_data.index[np.in1d(out_data.segment, segments)]:
        loc_idx = tuple(out_data.loc[i,['k','i', 'j']])
        out_data.loc[i,'flow_gw'] = str_recharge_data[loc_idx]
        out_data.loc[i,'flow_sw'] = str_flow_data[loc_idx]
        out_data.loc[i,'c_gw'] = gw_conc_data[loc_idx]

    #set inital surface water concentrations
    # discuss surface water concentration setting This is not non-complex
    # possibly set from the inital concentrations of MT3D (so as not to have to pass further arguments)
    for key in str_seg_trib_dict.keys():
        if str_seg_trib_dict[key] is not None:
            continue

        idx = out_data.index[(out_data.segment == key) & (out_data.reach ==1)]
        out_data.loc[idx,'c_sw'] = 0 # set influx of water to zero concentration #  a start for now fix with actual bc concentrations

    # calculate concentration data downstream
    seg_completed = pd.Series(index=range(1,45),data=False)
    temp_seg_list = range(1,45)
    while not all(seg_completed):
        for seg in temp_seg_list:
            if not (str_seg_trib_dict[seg] is None or all(seg_completed.loc[str_seg_trib_dict[seg]])):
                continue

            for re in range(1,out_data.reach[out_data.segment == seg].max()+1):
                if pd.notnull(out_data.c_sw.loc[(seg,re)]):
                    continue

                # special case for reach 1
                    #special case for overlapping #watch trib overlap behavior

                # special case for overlapping #watch trib overlap behavior

                # all others
                gwf = out_data.loc[(seg,re),'flow_gw']
                swf = out_data.loc[(seg,re),'flow_sw'] # this is flow out of the cell
                gwc = out_data.loc[(seg,re),'c_gw']
                swc_p = out_data.loc[(seg,re-1),'c_sw']

                if any(pd.isnull([gwf,swf,gwc,swc_p])):
                    raise ValueError ('null value for s/r: {}'.format((seg,re)))

                if gwf >= 0:  # if water is leaving stream to ground water then it won't change the stream concentration
                    out_data.loc[(seg,re),'c_sw'] = swc_p
                    continue

                #[] load mass balance
                sw_load = (swf-gwf)*swc_p + gwc*gwf

                out_data.loc[(seg,re,),'c_sw'] = sw_load/swf


                raise ValueError('not completed')
            seg_completed.loc[seg] = True

    return out_data


def get_all_segs(segment_id):
    # gets all tributary segments for a given feature
    segments = []
    segs_added = pd.Series(index=np.atleast_1d(segment_id),data=False)

    while not all(segs_added):
        for i in segs_added.index[np.invert(segs_added)]:
            segments.append(i)
            if str_seg_trib_dict[i] is not None:
                temp = pd.Series(index=str_seg_trib_dict[i],data=False)
                segs_added = segs_added.append(temp)
            segs_added.loc[i] = True

    return segments


if __name__ == '__main__':
    print(get_all_segs(18))