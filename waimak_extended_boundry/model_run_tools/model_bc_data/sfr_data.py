"""
 Author: Matt Hanson
 Created: 17/03/2020 8:26 AM
 """
from __future__ import division
import pandas as pd
import numpy as np
import os
import env

def get_base_sfr_data(): # todo check
    """
    get the pre-optimisation data, exlcude hcond etc????
    :return:
    """
    seg_data = pd.read_hdf(os.path.join(env.sdp_required,'sfr_data.hdf'),'seg_data')
    rch_data = pd.read_hdf(os.path.join(env.sdp_required,'sfr_data.hdf'),'reach_data')

    return seg_data, rch_data

def add_sfr_cond(model_id, rch_data): # todo
    """
    get teh sfr conductance for a given reach
    :param model_id:
    :return:
    """
    # todo possibly make netcdf of hcond valuse as there is some interpolation here... or see how long it takes...
    raise NotImplementedError

def add_sfr_inflows(model_id, seg_data): # todo
    """
    add inflow data to the segment data, the sfr inflows for the Cust, cust biwash and Eyre are parameterised while
    the inflows for the Waimakariri R. Ashley R. and Ashley Tribs, were held constant
    :param model_id:
    :return:
    """
    constant_inflows = {
        3: 11300000.0, # 'Waimakariri'
        5: 976320.0, # 'Ashley'
        8: 59213.0, # 'Glen Tui (Ashley)'
        15: 172535.0, # 'Garry (Ashley)
        22: 79632.0, # 'Bullock (Ashley)'
        24: 547213.0, # 'Okuku (Ashley)'
        28: 74527.0 # 'Makerikeri (Ashley)'
    }

    non_constant = {
        1: 171936.0, # 'Eyre'
        7: 14688.0, # 'Cust'
        19: 3456.0, # 'Race Bywash (Cust)'
    } # todo pull from paramsnetcdf

    raise NotImplementedError

def get_sfr_data(model_id, return_as_df=False):
    """
    get teh sfr data
    :param model_id: model id
    :param return_as_df: boolean if True redurnt as pd.DataFrame, otherwise return as record array for flopy
    :return:
    """
    seg_data, rch_data = get_base_sfr_data()

    seg_data = add_sfr_inflows(model_id, seg_data)
    rch_data = add_sfr_cond(model_id, rch_data)

    if return_as_df:
        return seg_data, rch_data
    else:
        seg_data = convert_seg_to_recarray(seg_data)
        rch_data = convert_rch_to_recarray(rch_data)
        return seg_data, rch_data


def convert_rch_to_recarray(rch_data):
    """
    convert teh reach data to teh record array needed for modflow
    :param rch_data:
    :return:
    """
    raise NotImplementedError

def convert_seg_to_recarray(seg_data):
    """
    convert the segment data to the record array needed for modflow
    :return:
    """
    raise NotImplementedError