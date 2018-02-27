"""
Author: matth
Date Created: 29/05/2017 11:52 AM
"""

from __future__ import division
import flopy
from drain_concentration import get_drn_samp_pts_dict, get_drn_samp_pts
import glob
import pandas as pd
import numpy as np
from get_str_rch_values import get_stream_duplication_array

# some changes in sfr require some changes here.

def streamflow_for_kskps(model_base_path, kstpkpers=None, drn_points=None, str_points=None):
    """
    get stream flows at points (specified by drn and str points for streams issued as drains and as stream respectivly)
    returns a pandas dataframe
    :param model_base_path: the base path to the model files (e.g. the name file without the extension)
    :param kstpkpers: either tuple (kstp, kper) or list of said tuples if none returns for all kstp,kper in the file
    :param drn_points: list of identifies for drain flow measurement points if None returns all avalible points) pass
                       an empty list to return none
    :param str_points: as above but for stream boundry conditions.
    :return: stream flows in model units pandas dataframe with multiindex of (kstp,kper) and columns of str and drn points
    """

    drn_samp_pts_dict = get_drn_samp_pts_dict()

    if drn_points is None:
        drn_points = get_drn_samp_pts()
    else:
        if not set(drn_points).issubset(get_drn_samp_pts()):
            raise ValueError('unexpected drn sample points: {}'.format(set(drn_points) - set(get_drn_samp_pts())))

    if str_points is None:
        str_points = get_str_sample_points()
    else:
        if not set(str_points).issubset(get_str_sample_points()):
            raise ValueError('unexpected str sample points: {}'.format(set(str_points) - set(get_str_sample_points())))

    drain_cbc = flopy.utils.CellBudgetFile('{}.cbd'.format(model_base_path))

    str_str_cbc = flopy.utils.CellBudgetFile('{}.csts'.format(model_base_path))

    str_model_cbc = flopy.utils.CellBudgetFile('{}.cstm'.format(model_base_path))

    if kstpkpers is None:
        kstpkpers = drain_cbc.get_kstpkper()

    kstpkpers = np.atleast_2d(kstpkpers)  # sort out the possiblity of only one kstpkper
    kstpkpers = [tuple(e) for e in kstpkpers]
    # set up multi index pandas dataframe (kstp, kper, columns=all the drain and stream points calculated)

    idx = pd.MultiIndex.from_tuples(kstpkpers, names=['kstp', 'kper'])
    outdata = pd.DataFrame(index=idx, columns=str_points + drn_points)
    for kstpkper in kstpkpers:
        kstp = kstpkper[0]
        kper = kstpkper[1]
        # create hydrograph data for drains
        drn_flow = drain_cbc.get_data(kstpkper=kstpkper, full3D=True)[0][0]
        for drn in drn_points:
            outdata.loc[kstp, kper].loc[drn] = _get_drn_flow(drn, drn_flow, drn_samp_pts_dict)

        # create hydrograph data for streams
        str_flow = str_str_cbc.get_data(kstpkper=kstpkper, full3D=True)[0][0]
        str_model_flow = str_model_cbc.get_data(kstpkper=kstpkper, full3D=True)[0][0]
        for stream in str_points:
            outdata.loc[kstp, kper].loc[stream] = _get_str_flow(stream, str_flow, str_model_flow, drn_flow,
                                                                drn_samp_pts_dict)
    return outdata


def _get_drn_flow(drn, drn_flows, drain_flow_dict):
    """
    function to get stream flow from drains or sd from stream to model flow
    :param drn:
    :param drn_flows:
    :param drain_flow_dict:
    :return:
    """
    mask_array = drain_flow_dict[drn]

    if drn_flows[mask_array].mask.sum() != 0:
        raise ValueError('masked values returned for {}'.format(drn))

    water_volume = drn_flows[mask_array].data.sum() * -1
    return water_volume


def _get_str_flow(stream, str_to_str_flow, str_to_model_flow, drn_flow, drn_samp_dict):
    """
    large wrapper to calculate any flows as stream points
    :param stream: the key for the location to calcuate the stream flow or sd
    :param str_to_str_flow: the 2d array of the stream to stream flow
    :param str_to_model_flow:  the 2d array of the stream to model flow
    :param drn_flow:  the 2d array of the drain to model flow
    :param drn_samp_dict: the drain sampling dictionary that links mask arrays to drain sample points
    :return:
    """

    str_function = _str_to_function_dict()[stream]
    water_volume = str_function(str_flow=str_to_str_flow, drn_flow=drn_flow, drn_samp_dict=drn_samp_dict,
                                str_model_flow=str_to_model_flow)

    return water_volume


def get_str_sample_points():
    """
    get the keys for the set up stream functions
    :return:
    """
    temp = _str_to_function_dict()
    return temp.keys()


def _str_to_function_dict():
    """
    map stream functions (below) to the stream name
    :return:
    """
    str_to_function = {
        'eyre_warren_hyg': _calc_eyre_warren_hydrograph,
        'ashley_sh1_hyg': _calc_ashley_sh1_hydrograph,
        'cust_therelkelds_hyg': _calc_cust_therelkelds_hydrograph,
        'cust_oxfordrd_hyg': _calc_cust_oxfordrd_hydrograph,
        'num7drn_cust_hyg': _calc_num7drain_hydrograph,
        'kaiapoi_waimak_hyg': _calc_kaiapoi_waimak_hydrograph,
        'eyre_swaz_sd': _calc_eyre_swaz_sd,
        'ashley_swaz_sd': _calc_ashley_swaz_sd,
        'custmain_swaz_sd': _calc_custmain_swaz_sd,
        'custr_swaz_sd': _calc_custr_swaz_sd,
        'num7drn_swaz_sd': _calc_num7drn_swaz_sd
    }
    return str_to_function



def _get_str_flow_at_pt(str_flow, row, col):
    """
    internal function to deal with the challenges of duplication of stream points and return the correct cell by cell
    flow works on either the stream flow or the model to stream flow
    :param str_flow:
    :param row:
    :param col:
    :return:
    """
    dup = get_stream_duplication_array()[row, col]
    if dup == 0:
        raise ValueError('duplication array returned zero at row: {} col:{}'.format(row, col))
    water_vol, mask = str_flow.data[row, col], str_flow.mask[row, col]
    if mask:
        raise ValueError('retuned masked values')
    return water_vol / dup


# functions to calcualte the stream flow or stream depletion at each of a number of sites

def _calc_eyre_warren_hydrograph(str_flow, **kwargs):
    """
    This is the first of many stream functions to calculate the flow or stream depletion at a given point.
    :param str_flow:
    :param kwargs: not used, to prevent errors when passing arguments across all stream functions
    :return:
    """
    water_vol = _get_str_flow_at_pt(str_flow, 56, 119)
    return water_vol


def _calc_ashley_sh1_hydrograph(str_flow, **kwargs):
    water_vol = _get_str_flow_at_pt(str_flow, 32, 309)
    return water_vol


def _calc_cust_therelkelds_hydrograph(str_flow, drn_flow, drn_samp_dict, **kwargs):
    water_vol = _get_str_flow_at_pt(str_flow, 80, 281)
    n7drn_flow = _get_drn_flow('num7drain_swaz', drn_flow, drn_samp_dict)
    return water_vol + n7drn_flow


def _calc_cust_oxfordrd_hydrograph(str_flow, **kwargs):
    water_vol = _get_str_flow_at_pt(str_flow, 53, 220)

    return water_vol


def _calc_num7drain_hydrograph(str_flow, drn_flow, drn_samp_dict, **kwargs):
    water_vol = _get_str_flow_at_pt(str_flow, 77, 280)

    n7drn_flow = _get_drn_flow('num7drain_swaz', drn_flow, drn_samp_dict)
    return water_vol + n7drn_flow


def _calc_kaiapoi_waimak_hydrograph(str_flow, drn_flow, drn_samp_dict, **kwargs):
    water_vol = _get_str_flow_at_pt(str_flow, 86, 290)  # at bottom of cust stream bcs (exclueds part of num7drn)
    drns_to_add = ['num7drain_swaz', 'cam_kaiapoi', 'courtenay_kaiapoi', 'kaiapoi_swaz', 'ohoka_swaz']
    for drn in drns_to_add:
        water_vol += _get_drn_flow(drn, drn_flow, drn_samp_dict)

    return water_vol

# by round 2 fix the overlap problem.

def _calc_eyre_swaz_sd(str_model_flow, drn_samp_dict, **kwargs):
    """
    the first of the stream depltion calculations
    :param str_model_flow:
    :param drn_samp_dict:
    :param kwargs:
    :return:
    """
    water_vol = _get_drn_flow('str_eyre_swaz', str_model_flow, drn_samp_dict)

    return water_vol


def _calc_ashley_swaz_sd(str_model_flow, drn_samp_dict, drn_flow, **kwargs):
    water_vol = _get_drn_flow('str_ashley_swaz', str_model_flow, drn_samp_dict)
    water_vol += _get_drn_flow('ashley_swaz', drn_flow, drn_samp_dict)

    return water_vol


def _calc_custmain_swaz_sd(str_model_flow, drn_samp_dict, **kwargs):
    water_vol = _get_drn_flow('str_custmaindrn_swaz', str_model_flow, drn_samp_dict)

    return water_vol


def _calc_custr_swaz_sd(str_model_flow, drn_samp_dict, **kwargs):
    water_vol = _get_drn_flow('str_cust_swaz', str_model_flow, drn_samp_dict)

    return water_vol


def _calc_num7drn_swaz_sd(str_model_flow, drn_samp_dict, drn_flow, **kwargs):
    water_vol = _get_drn_flow('str_num7drain_swaz', str_model_flow, drn_samp_dict)
    water_vol += _get_drn_flow('num7drain_swaz', drn_flow, drn_samp_dict)

    return water_vol


# the below are depreciated
def get_str_targets():
    # depreciated but held for past compatability
    str_targets = pd.DataFrame(
        {'name': ['cust_skewbridge', 'cust_threlkelds', 'ashley_sh1', 'eyre_waimak', 'cust_kaiapoi', 'eyre_wolffs'],
         'seg': [37, 36, 34, 17, 37, np.nan],
         'reach': [12, 4, 15, 53, 14, np.nan],
         'row': [85, 81, 32, 120, 86, 67],
         'col': [289, 282, 309, 256, 290, 152],
         'overlap': [0, 1, 0, 0, 0, np.nan]})
    str_targets = str_targets.set_index('name')
    return str_targets


def create_hydrographs_for_transient(model_dir):
    """
    function to look at the transient runs, kind of hackish sorry future me depreciated, but held for past compataiblity
    :param model_dir:
    :return:
    """
    drn_samp_pts_dict = get_drn_samp_pts_dict()
    str_vals = get_str_targets()
    drn_points = [
        'cam_youngs',
        'taranaki_gressons',
        'waikuku_sh1',
        'taranaki_preeces',
        'saltwater_factory',
        'kaiapoi_island',
        'ohoka_island',
        'kaiapoi_harpers',
    ]

    stress_to_month = {0: None,
                       1: 7,
                       2: 8,
                       3: 9,
                       4: 10,
                       5: 11,
                       6: 12,
                       7: 1,
                       8: 2,
                       9: 3,
                       10: 4,
                       11: 5,
                       12: 6}

    str_points = ['cust_skewbridge', 'cust_threlkelds', 'ashley_sh1']

    drn_cvc_path = glob.glob('{}/*.cbd'.format(model_dir))
    if len(drn_cvc_path) > 1:
        raise ValueError('more paths than expected in {}'.format(model_dir))
    drn_cvc_path = drn_cvc_path[0]

    str_cvc_path = glob.glob('{}/*.csts'.format(model_dir))
    if len(str_cvc_path) > 1:
        raise ValueError('more paths than expected in {}'.format(model_dir))
    str_cvc_path = str_cvc_path[0]

    outdata = pd.DataFrame(index=range(1, 13), columns=drn_points + str_points)
    drn_cell_flow = flopy.utils.CellBudgetFile(drn_cvc_path)
    str_cell_flow = flopy.utils.CellBudgetFile(str_cvc_path)
    # create hydrograph data for drains
    # average the two values in the month as a first stab
    for sper in range(1, 13):
        month = stress_to_month[sper]
        first_step_drn = drn_cell_flow.get_data(kstpkper=(0, sper), full3D=True)[0][0]
        second_step_drn = drn_cell_flow.get_data(kstpkper=(1, sper), full3D=True)[0][0]
        for drn in drn_points:
            mask_array = drn_samp_pts_dict[drn]

            if first_step_drn[mask_array].mask.sum() != 0:
                raise ValueError('masked values returned for {}'.format(drn))

            water_volume1 = first_step_drn[mask_array].data.sum()
            water_volume2 = second_step_drn[mask_array].data.sum()
            outdata.loc[month, drn] = (water_volume1 + water_volume2) / -2

        # create hydrograph data for streams
        first_str_step = str_cell_flow.get_data(kstpkper=(0, sper), full3D=True)[0][0]
        second_str_step = str_cell_flow.get_data(kstpkper=(1, sper), full3D=True)[0][0]
        for stream in str_points:
            row, col, modval = str_vals.loc[stream, ['row', 'col', 'overlap']]
            modval += 1
            water_volume1 = first_str_step[row, col]
            water_volume2 = second_str_step[row, col]
            outdata.loc[month, stream] = (water_volume1 + water_volume2) / 2
    return outdata


if __name__ == '__main__':
    test = streamflow_for_kskps(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\base_model_runs\base_stream_depletion_mf\m_strong_vert-fully_naturalized_sd_run\m_strong_vert-fully_naturalized_sd_run")
    test.to_csv(r"C:\Users\MattH\Downloads\test_str_flow.csv")
    print('done')
