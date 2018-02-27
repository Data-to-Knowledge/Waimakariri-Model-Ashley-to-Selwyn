"""
Author: matth
Date Created: 23/05/2017 1:48 PM
"""

from __future__ import division
from core import env
from m_wraps.base_modflow_wrapper import get_base_mf_ss
import pickle
import os
from users.MH.Waimak_modeling.supporting_data_path import sdp
from polygon_to_model_array import shape_file_to_model_array
import numpy as np
import pandas as pd
from basic_tools import df_to_array
from drain_concentration import get_drn_samp_pts_dict

# these are only for m_strong_vert
def get_base_str(recalc=False):
    picklepath = '{}/inputs/pickeled_files/base_stream_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        base_str = pickle.load(open(picklepath))
        return base_str

    org_m = get_base_mf_ss()
    base_str = org_m.str.stress_period_data.data[0]
    pickle.dump(base_str, open(picklepath, 'w'))

    return base_str

def get_ibound(recalc=False):
    picklepath = '{}/inputs/pickeled_files/ibound.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        ibound = pickle.load(open(picklepath))
        return ibound

    org_m = get_base_mf_ss()
    ibound = org_m.bas6.ibound.array
    pickle.dump(ibound, open(picklepath, 'w'))

    return ibound


def get_base_seg_data(recalc=False):
    picklepath = '{}/inputs/pickeled_files/base_stream_segments.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        base_str = pickle.load(open(picklepath))
        return base_str

    org_m = get_base_mf_ss()
    base_str = org_m.str.segment_data[0]
    pickle.dump(base_str, open(picklepath, 'w'))

    return base_str


def get_base_rch(recalc=False):
    picklepath = '{}/inputs/pickeled_files/base_recharge_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        base_rch = pickle.load(open(picklepath))
        return base_rch

    org_m = get_base_mf_ss()
    base_rch = org_m.rch.rech.array[0, 0]
    pickle.dump(base_rch, open(picklepath, 'w'))

    return base_rch


def get_gmp_rch(recalc=False):
    picklepath = '{}/inputs/pickeled_files/gmp_recharge_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        gmp_rch = pickle.load(open(picklepath))
        return gmp_rch

    base_rch = get_base_rch()

    factor = shape_file_to_model_array("{}/inputs/shp_files/rch_s/gmp_to_cmp.shp".format(sdp),
                                       'OMPLSRFact', alltouched=True)
    factor[np.isnan(factor)] = 1
    factor[~np.isclose(factor,
                       1)] *= 1.3  # this increase is to deal with the differences between fouad's Cmp layer and all other cmp layers

    gmp_rch = base_rch / factor
    pickle.dump(gmp_rch, open(picklepath, 'w'))

    return gmp_rch


def get_true_cmp_rch(recalc=False):
    picklepath = '{}/inputs/pickeled_files/true_cmp_recharge_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        true_cmp_rch = pickle.load(open(picklepath))
        return true_cmp_rch

    # the current model has something that is probably closer to gmp than cmp, so here we back calculate a more
    # realistic cmp layer.  the naming true CMP is to not kill my naming.
    base_rch = get_base_rch()

    factor = shape_file_to_model_array("{}/inputs/shp_files/rch_s/gmp_to_cmp.shp".format(sdp),
                                       'OMPLSRFact', alltouched=True)
    factor[np.isnan(factor)] = 1

    true_cmp_rch = base_rch * factor
    pickle.dump(true_cmp_rch, open(picklepath, 'w'))
    return true_cmp_rch


def get_nat_rch(recalc=False):
    picklepath = '{}/inputs/pickeled_files/nat_recharge_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        nat_rch = pickle.load(open(picklepath))
        return nat_rch

    base_rch = get_base_rch()
    nat_rch = shape_file_to_model_array("{}/inputs/shp_files/rch_s/drylandlsr_nodes_0.shp".format(sdp),
                                        'Topography', alltouched=True)

    nat_rch[nat_rch <= 0] = np.nan
    nat_rch *= 2.73973e-6
    idx = np.isnan(nat_rch)
    nat_rch[idx] = base_rch[idx]
    pickle.dump(nat_rch, open(picklepath, 'w'))
    return nat_rch


def get_stream_duplication_array(recalc=False):
    # it looks like things are not working perfectly with this... really tease this out for N
    picklepath = '{}/inputs/pickeled_files/str_dup_array.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        str_dup_array = pickle.load(open(picklepath))
        return str_dup_array

    base_str = pd.DataFrame(get_base_str())
    str_dup_array = np.zeros((190, 365))
    for i in base_str.index:
        row, col = base_str.loc[i, ['i', 'j']]
        str_dup_array[row, col] += 1

    pickle.dump(str_dup_array, open(picklepath, 'w'))
    return str_dup_array


def get_stream_seg_dict(recalc=False):
    picklepath = '{}/inputs/pickeled_files/str_seg_dict.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        str_seg_dict = pickle.load(open(picklepath))
        return str_seg_dict

    base_str = pd.DataFrame(get_base_str())

    keys = ['str_ashley_swaz',
            'str_cust_swaz', 'str_custmaindrn_swaz', 'str_eyre_swaz',
            'str_num7drain_swaz']
    str_seg_dict = {}
    seg_array = df_to_array(base_str, 'segment')
    drn_dict = get_drn_samp_pts_dict()
    temp = []
    for key in keys:
        t = set(seg_array[drn_dict[key]])
        str_seg_dict[key] = t
        temp.extend(t)
    str_seg_dict['str_waimak_swaz'] = set(base_str['segment']) - set(temp)
    pickle.dump(str_seg_dict, open(picklepath, 'w'))
    return str_seg_dict


def aqualinc_seg_dict():
    aslink = {
        'a1': [5, 8],
        'a2': [12, 15, 18, 22, 30],
        'a3': [24, 31],
        'a4': [28, 32],
        'a5': [33],
        'a6': [34],
        'c1': [7, 10],
        'c2': [13],
        'c3': [16],
        'c4': [19, 21],
        'c5': [23, 42, 43, 26, 44, 36, 37],
        'e1': [1],
        'e3': [2],
        'e4': [4],
        'e5': [11],
        'e6': [17],
        'w1': [3],
        'w2': [6],
        'w3': [9],
        'w4': [14],
        'w5': [20, 29],
        'w6': [35],
    }
    return aslink


def get_base_drn_cells(recalc=False):
    picklepath = '{}/inputs/pickeled_files/base_drn_spd.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        base_drn = pickle.load(open(picklepath))
        return base_drn

    org_m = get_base_mf_ss()
    base_drn = org_m.drn.stress_period_data.data[0]
    pickle.dump(base_drn, open(picklepath, 'w'))

    return base_drn


if __name__ == '__main__':
    get_stream_seg_dict(True)
