"""
Author: matth
Date Created: 14/06/2017 12:27 PM
"""

from __future__ import division
import numpy as np
import pandas as pd
import flopy_mh as flopy
from users.MH.Waimak_modeling.model_tools.get_str_rch_values import get_base_str, get_base_seg_data, get_stream_seg_dict, aqualinc_seg_dict, get_base_drn_cells
from users.MH.Waimak_modeling.supporting_data_path import sdp
import geopandas as gpd
from copy import deepcopy
import pickle
import os
from waimak_extended_boundry import smt
from stream_elvs import get_reach_elv


def create_sfr_package(m, version=1, seg_v=1, reach_v=1):
    """
    wrapper to create sfr packages
    :param m: the model on which to create the package
    :param version: which sfr version to use.  implemented here to record multiple attempts at the sfr
    :param seg_v: which version of the segment data to use here to record multiple attempts at the sfr
    :param reach_v: which verson of the reach data to use here to record multiple attempts at the sfr
    :return:
    """
    if version == 1:
        _create_sfr_version_1(m, seg_v, reach_v)


def _create_sfr_version_1(m, seg_v, reach_v):
    """
    sfr version one, no transient flow routing, no unsaturated zone flow
    :param m: modflow model object on which to build the package
    :return:
    """
    reach_data = _get_reach_data(reach_v)
    segment_data = _get_segment_data(seg_v)
    sfr = flopy.modflow.ModflowSfr2(m,
                                    nstrm=len(reach_data),
                                    nss=len(segment_data),
                                    nsfrpar=0,  # this is not implemented in flopy not using
                                    nparseg=0,  # this is not implemented in flopy not using
                                    const=86400,
                                    dleak=0.0001,  # default, likely ok
                                    ipakcb=740,
                                    istcb2=0,
                                    isfropt=1,  # Not using as nstrm >0 therfore no Unsaturated zone flow
                                    nstrail=10,  # not using as no unsaturated zone flow
                                    isuzn=1,  # not using as no unsaturated zone flow
                                    nsfrsets=30,  # not using as no unsaturated zone flow
                                    irtflg=0,  # no transient sf routing
                                    numtim=2,  # Not using transient SFR
                                    weight=0.75,  # not using if irtflg = 0 (no transient SFR)
                                    flwtol=0.0001,  # not using if irtflg = 0 (no transient SFR)
                                    reach_data=reach_data,
                                    segment_data=segment_data,
                                    channel_geometry_data=None,  # using rectangular profile so not used
                                    channel_flow_data=None,  # icalc = 1 or 2 so not using
                                    dataset_5={0:[len(segment_data),0,0]},
                                    reachinput=True,  # using for reach input
                                    transroute=False,  # no transient sf routing
                                    tabfiles=False,  # not using
                                    tabfiles_dict=None,  # not using
                                    unit_number=717)


def _get_segment_data(seg_v):
    """
    wrapper to get the segment data
    :param seg_v: the version here to record multiple attempts at the sfr
    :return:
    """
    if seg_v == 1:
        seg_data = _seg_data_v1()
    else:
        raise ValueError('unexpected version for segment: {}'.format(seg_v))
    return seg_data


def _get_reach_data(reach_v):
    """
    wrapper to get the reach data
    :param reach_v: the version to use here to record multiple attempts at the sfr
    :return:
    """
    if reach_v == 1:
        reach_data = _reach_data_v1()
    else:
        raise ValueError('unexpected version for reach: {}'.format(reach_v))
    return reach_data


def _reach_data_v1(recalc=False):
    """
    load reach version one (the only version used the extended boundary waimakariri proces as of 20/10/2017)
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return: modflow sfr stress period data.  numpy record array
    """
    pickle_path = '{}/sfr_reach_v1.p'.format(smt.pickle_dir)

    if os.path.exists(pickle_path) and not recalc:
        reach_data = pickle.load(open(pickle_path))
        return reach_data

    temp_str_data = _get_base_stream_values()

    outdata = flopy.modflow.ModflowSfr2.get_empty_reach_data(len(temp_str_data.index),default_value=0)
    data_to_pass = {'i': 'i', 'j': 'j', 'k': 'k', 'reach': 'ireach', 'segment': 'iseg', 'slope': 'slope',
                    'stop': 'strtop'}

    # fill drain and link cells
    temp_str_data.loc[pd.isnull(temp_str_data['slope']),'slope'] = temp_str_data['slope'][temp_str_data['segment']==27].mean()
    temp_str_data.loc[pd.isnull(temp_str_data['stop']),'stop'] = temp_str_data.loc[pd.isnull(temp_str_data['stop']),'elev']
    temp_str_data.loc[temp_str_data['segment'] == 27,'stop'] = temp_str_data.loc[temp_str_data['segment']==27,'stop'].interpolate()

    for key in data_to_pass.keys():
        outdata[data_to_pass[key]] = temp_str_data[key]

    outdata = _define_reach_length(outdata, mode='cornering')
    outdata['strthick'] = 1

    # define k from the conductance length and width as a starting point
    seg_data = pd.DataFrame(_seg_data_v1())
    temp_str_data.loc[temp_str_data['segment'] == 27, 'cond'] = 500  # just set as a default value for now
    temp_str_data['width'] = np.nan
    for seg, width in zip(seg_data['nseg'],seg_data['width1']):
        temp_str_data.loc[temp_str_data['segment']==seg, 'width'] = width
    outdata['strhc1'] = temp_str_data.loc[:,'cond']/(outdata['rchlen'] * temp_str_data['width'])

    outdata['strtop'] = np.array(get_reach_elv()) - 0.5

    # fix wierd segment numbering for reaches
    temp_str_data = pd.DataFrame(outdata)
    old_to_temp_seg = dict([[38,125],[42,127],[41,126],[39,136],[43,137],[25,138],[26,139],[44,141],[27,142],[36,143],[37,144]])
    temp_str_data = temp_str_data.replace({'iseg':old_to_temp_seg})
    temp_to_new_seg = dict([[125,25],[127,27],[126,26],[136,36],[137,37],[138,38],[139,39],[141,41],[142,42],[143,43],[144,44]])
    temp_str_data = temp_str_data.replace({'iseg':temp_to_new_seg})

    #fix one backwards flowing reach
    temp_str_data.loc[temp_str_data.iseg==19,'strtop'] = 72.44

    outdata = temp_str_data.to_records(False).astype(flopy.modflow.ModflowSfr2.get_default_reach_dtype())


    elv = smt.calc_elv_db()
    temp = pd.DataFrame(outdata)
    str_tops = smt.df_to_array(temp,'strtop')
    if any((str_tops > elv[0]).flatten()):
        raise ValueError('streams with elevation above surface')

    if any((str_tops-1 <= elv[1]).flatten()):
        raise ValueError('streams below layer 1')


    pickle.dump(outdata,open(pickle_path,'w'))
    return outdata

def _seg_data_v1(recalc=False):
    """
    get segment version 1 the only version used in the extended waimakariri model as of 20/10/2017
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return: sfr seg data (np.record array)
    """
    pickle_path = '{}/sfr_seg_v1.p'.format(smt.pickle_dir)

    if os.path.exists(pickle_path) and not recalc:
        seg_data = pickle.load(open(pickle_path))
        return seg_data

    seg_data = flopy.modflow.ModflowSfr2.get_empty_segment_data(44, default_value=0)
    seg_data['nseg'] = range(1, 45)
    seg_data['icalc'] = 1

    str_data = _get_base_stream_values()
    str_data.loc[pd.isnull(str_data['flow']),'flow'] = 0
    str_data.loc[str_data['flow']<0,'flow'] = 0

    # outseg
    org_seg_data = pd.DataFrame(get_base_seg_data())
    org_seg_data['seg'] = range(1,len(org_seg_data)+1)
    trib_dict = dict(np.array(org_seg_data[['itrib01','seg']]))
    trib_dict.update(np.array(org_seg_data[['itrib02','seg']]))
    trib_dict.update(np.array(org_seg_data[['itrib03','seg']]))
    for i,seg in enumerate(seg_data['nseg']):
        if seg in trib_dict.keys():
            seg_data['outseg'][i] = trib_dict[seg]
    # get flows
    for i,seg in enumerate(seg_data['nseg']):
        seg_data['flow'][i] = str_data['flow'][(str_data['reach']==1) & (str_data['segment']==seg)].iloc[0]
        # take the max roughness of segment 11 (lower eyre river) all others only have 1 roughness
        seg_data['roughch'][i] = str_data['rough'][str_data['segment'] == seg].max()


    # hard variable
    seg_dict = get_stream_seg_dict()

    # the aqualinc report suggests that the eyre and ashley are relativly constant width so set all width to the value
    # provided # this did not work for some reason... and instead teh widths were set to 1
    seg_data['width1'][np.in1d(seg_data['nseg'],seg_dict['str_ashley_swaz'])] = 20  # set ashley width to 20 m from aqualinc
    seg_data['width2'][np.in1d(seg_data['nseg'],seg_dict['str_ashley_swaz'])] = 20  # set ashley width to 20 m from aqualinc
    seg_data['width1'][np.in1d(seg_data['nseg'],seg_dict['str_eyre_swaz'])] = 5.5  # set eyre width to 5.5 m from aqualinc
    seg_data['width2'][np.in1d(seg_data['nseg'],seg_dict['str_eyre_swaz'])] = 5.5  # set eyre width to 5.5 m from aqualinc


    aq_data = pd.read_excel("{}\inputs\str_data_from_aqualinc.xlsx".format(sdp))
    aq_data = aq_data.set_index('Vaue')
    seg_aq_dict = aqualinc_seg_dict()
    # the aqualinc report suggests that the cust and the waimak vary substantially across their length, so I have set
    # them accordingly
    for key in seg_aq_dict.keys():
        if 'a' in key or 'e' in key: # skip the ashley and the eyre
            continue
        segs =  seg_aq_dict[key]
        wdth = aq_data.loc['width',key]
        seg_data['width1'][np.in1d(seg_data['nseg'],segs)] = wdth
        seg_data['width2'][np.in1d(seg_data['nseg'],segs)] = wdth

    # the cust tributaries are not detailed in the aqualinc report, I'll define these by approximation from areal imagry
    seg_data['width1'][np.isclose(seg_data['width1'], 0)] = 1 # define all undefined widths as 1 m these should just be the drains
    seg_data['width2'][np.isclose(seg_data['width2'], 0)] = 1 #define all undefined widths as 1 m

    #fix inconsistant ordering
    seg_data = pd.DataFrame(seg_data)
    old_to_temp_seg = dict([[38,125],[42,127],[41,126],[39,136],[43,137],[25,138],[26,139],[44,141],[27,142],[36,143],[37,144]])
    seg_data = seg_data.replace({'nseg':old_to_temp_seg,'outseg':old_to_temp_seg})
    temp_to_new_seg = dict([[125,25],[127,27],[126,26],[136,36],[137,37],[138,38],[139,39],[141,41],[142,42],[143,43],[144,44]])
    seg_data = seg_data.replace({'nseg':temp_to_new_seg,'outseg':temp_to_new_seg}).sort_values('nseg')
    seg_data.loc[18,'flow'] = 3456
    seg_data = seg_data.to_records(False).astype(flopy.modflow.ModflowSfr2.get_default_segment_dtype())

    pickle.dump(seg_data,open(pickle_path,'w'))
    return seg_data


def _get_base_stream_values():
    """
    get all current data and locations for the stream segments.  it will be a mix of drains and stream and made up data
    :return:
    """
    # get rid of duplicates from the orginal stream package
    str_data = pd.DataFrame(get_base_str())
    dup = str_data.duplicated(['i', 'j'], False)
    str_data = str_data[(~dup) | (str_data['reach'] == 1)]

    # add drain cells and additional cells that we need
    temp_drn_data = gpd.read_file("{}/inputs/shp_files/num7_drn_to_str/drn_to_str.shp".format(sdp))
    temp_drn_data = temp_drn_data[['elv', 'lat', 'lon', 'reach', 'xdata', 'ydata']]
    temp_drn_data = temp_drn_data.rename(columns={'xdata': 'j', 'ydata': 'i'})
    org_drn = pd.DataFrame(get_base_drn_cells())
    temp = pd.concat((org_drn, temp_drn_data))
    drn_data = temp[temp.duplicated(['i','j'],keep='last')] # transition to the actual drain data rather than the GIS data, which is not what i need
    drn_data.loc[:,'reach'] = temp_drn_data.loc[:,'reach'].values

    # add other link cells
    idx = len(drn_data.index)
    for i, row in enumerate(range(64, 69)):
        drn_data.loc[idx + i, ['i', 'j', 'reach']] = (row, 277, idx + i + 1)
    drn_data['segment'] = 27
    drn_data['k'] = 0

    # combine drn cells and str_data
    str_data['reach'][str_data['segment'] == 27] += drn_data['reach'].max()
    str_data = pd.concat((drn_data, str_data))
    str_data = str_data.sort_values(['segment', 'reach']).reset_index()

    return str_data


def _define_reach_length(reach_data, mode='cornering'):
    """
    define the reach length for the SFR package
    :param reach_data: dataframe of reach data
    :param mode: the mode to use to determine reach length options:
                 constant: use the length of the grid cell (e.g. 200 m)
                 cornering: use the length of the grid cell for straight segments and use hypotenuse (of the half cell)
                            for corner cells assume kitty corner to kitty corner cells move straight across the cell
                            (e.g. the hypotius of the full cell) for kitty corner to adjacent assume that the stream
                            goes to teh center of the cell (along hypotenus of the half cell) and then straigth out to
                            the next cell (half cell lenght)

    :return: rech_data with length filled in
    """
    wrd = deepcopy(reach_data)
    cell_dim = 200
    if mode == 'constant':
        wrd['rchlen'] = cell_dim
    elif mode == 'cornering':
        rchlen = np.zeros(len(wrd['rchlen']))
        for i, rch, seg in zip(range(0, len(rchlen)), wrd['ireach'], wrd['iseg']):
            if rch == 1:  # assume full length for the first segment of each reach
                rchlen[i] = cell_dim
                continue
            if rch == wrd['ireach'][wrd['iseg'] == seg].max():  # assume full length for the last reach of the segment
                rchlen[i] = cell_dim
                continue

            # for all middle reaches check to see if it is a corner and then assign values associated with that
            c_row, c_col, c_seg, c_rch = wrd[['j', 'i', 'iseg', 'ireach']][i]
            p_row, p_col, p_seg, p_rch = wrd[['j', 'i', 'iseg', 'ireach']][i - 1]
            n_row, n_col, n_seg, n_rch = wrd[['j', 'i', 'iseg', 'ireach']][i + 1]

            # some checks to make sure everything is working right
            if p_seg != seg or n_seg != seg:
                raise ValueError('different segments adjacent for item: {}, reach: {}, seg: {} '.format(i, rch, seg))
            elif n_rch != rch + 1:
                raise ValueError('unexpected reaches adjacent for item: {}, reach: {}, seg: {} '.format(i, rch, seg))
            elif p_rch != rch - 1:
                raise ValueError('unexpected reaches adjacent for item: {}, reach: {}, seg: {} '.format(i, rch, seg))
            elif not np.in1d([p_row,n_row], range(c_row-1,c_row+2)).all():
                raise ValueError('rows not adjacent for item: {}, reach: {}, seg: {} '.format(i, rch, seg))
            elif not np.in1d([p_col,n_col], range(c_col-1,c_col+2)).all():
                raise ValueError('cols not adjacent for item: {}, reach: {}, seg: {} '.format(i, rch, seg))

            # different options for based on the location of the previous and next reach
            if p_row==c_row==n_row:  # straight segment across rows
                rchlen[i] = cell_dim
            elif p_col==c_col==n_col:  # straight segment across columns
                rchlen[i] = cell_dim
            elif p_col == c_col:
                if c_row == n_row: # assume the river cuts the corner
                    rchlen[i] = ((cell_dim/2)**2 + (cell_dim/2)**2)**0.5
                else: # assume the river goes to the centre of the cell and out
                    rchlen[i] = ((cell_dim)**2 + (cell_dim)**2) ** 0.5 + cell_dim/2
            elif p_row == c_row:
                if c_col == n_col: # assume the river cuts the corner
                    rchlen[i] = ((cell_dim / 2)**2 + (cell_dim / 2)**2) ** 0.5
                else: # assume the river goes to the centre of the cell and out
                    rchlen[i] = ((cell_dim)**2 + (cell_dim)**2) ** 0.5 + cell_dim/2
            else:
                if c_row == n_row or c_col == n_col: # assume the river goes to the centre of the cell and out
                    rchlen[i] = ((cell_dim / 2)**2 + (cell_dim / 2)**2) ** 0.5
                else: # assume the river cuts acorss the entire cell linearly
                    rchlen[i] = ((cell_dim)**2 + (cell_dim)**2) ** 0.5 + cell_dim/2
        wrd['rchlen'] = rchlen
    else:
        raise ValueError('unexpected argurment for mode: {}'.format(mode))
    return wrd

if __name__ == '__main__':
    # tests
    save = True
    seg = pd.DataFrame(_seg_data_v1(False))
    reach = pd.DataFrame(_reach_data_v1(False))
    if save:
        reach = smt.add_mxmy_to_df(reach)
        # this has been switched around reach.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\raw_sw_samp_points\sfr\all_sfr.csv")
        seg.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\raw_sw_samp_points\sfr\seg_data.csv")

    else:
        elv = smt.calc_elv_db()
        g=reach.groupby(['iseg'])
        bot = g.aggregate({'strtop':np.min})
        top = g.aggregate({'strtop':np.max})
        problems = []
        for i in seg.index:
            segment, outseg = seg.loc[i,['nseg','outseg']]
            if outseg ==0:
                continue
            ttop = top.loc[outseg, 'strtop']
            tbot = bot.loc[segment, 'strtop']
            if ttop>tbot:
                problems.append((segment,outseg))

        print(problems)
