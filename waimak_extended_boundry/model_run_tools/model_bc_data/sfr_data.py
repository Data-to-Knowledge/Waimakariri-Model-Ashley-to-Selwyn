"""
 Author: Matt Hanson
 Created: 17/03/2020 8:26 AM
 """
from __future__ import division
import pandas as pd
import numpy as np
import os
import env
import netCDF4 as nc


def get_base_sfr_data():  # todo check
    """
    get the pre-optimisation data, exlcude hcond etc????
    :return:
    """
    seg_data = pd.read_hdf(os.path.join(env.sdp_required, 'sfr_data.hdf'), 'seg_data')
    rch_data = pd.read_hdf(os.path.join(env.sdp_required, 'sfr_data.hdf'), 'reach_data')

    return seg_data, rch_data


def add_sfr_cond(model_id, rch_data):  # todo check
    """
    get teh sfr conductance for a given reach
    :param model_id:
    :return:
    """

    rch_data = rch_data.copy()
    param_data = nc.Dataset(os.path.join(env.sdp_required, "nsmc_params_obs_metadata.nc"))

    nsmc_num = int(model_id[-6:])
    nidx = np.where(param_data.variables['nsmc_num'][:] == nsmc_num)[0][0]
    conduct_names = param_data.variables['sfr_cond'][:]
    conduct_vals = param_data.variables['sfr_cond_val'][nidx][:]
    cond_dict = {k:v for k,v in zip(conduct_names,conduct_vals)}

    # get the parameter mapper
    param_mapper = pd.read_table(os.path.join(env.sdp_required,'base_for_nsmc_real/sfr_segdata.tpl'),
                                 skiprows=1, index_col=0)

    for nseg, hcond1_nm, hcond2_nm in param_mapper[['hcond1','hcond2']].itertuples(True,None):
        hcond1 = cond_dict[hcond1_nm]
        hcond2 = cond_dict[hcond2_nm]
        m = (hcond2 - hcond1) / rch_data.loc[rch_data.iseg==nseg,'rchlen'].sum()  # slope of conductance down segment
        # Series of conductance linearly varying from upper to lower
        rhcond = (m * (rch_data.loc[rch_data.iseg==nseg,'rchlen'].cumsum() -
                       rch_data.loc[rch_data.iseg==nseg,'rchlen'].rchlen / 2.0)) \
                 + hcond1  # m*x+c - x is distance along segment to cell midpoit

        rch_data.loc[rch_data.iseg == nseg, 'strhc1'] = rhcond

    return rch_data


def add_sfr_inflows(model_id, seg_data):  # todo check
    """
    add inflow data to the segment data, the sfr inflows for the Cust, cust biwash and Eyre are parameterised while
    the inflows for the Waimakariri R. Ashley R. and Ashley Tribs, were held constant
    :param model_id:
    :return:
    """
    seg_data = seg_data.copy()
    param_data = nc.Dataset(os.path.join(env.sdp_required, "nsmc_params_obs_metadata.nc"))

    nsmc_num = int(model_id[-6:])
    nidx = np.where(param_data.variables['nsmc_num'][:] == nsmc_num)[0][0]

    eyre = float(param_data.variables['top_e_flo'][nidx])
    cust = float(param_data.variables['top_c_flo'][nidx])
    race_biwash = float(param_data.variables['mid_c_flo'])

    inflows = {
        3: 11300000.0,  # 'Waimakariri'
        5: 976320.0,  # 'Ashley'
        8: 59213.0,  # 'Glen Tui (Ashley)'
        15: 172535.0,  # 'Garry (Ashley)
        22: 79632.0,  # 'Bullock (Ashley)'
        24: 547213.0,  # 'Okuku (Ashley)'
        28: 74527.0,  # 'Makerikeri (Ashley)'

        # non constant
        1: eyre,  # 'Eyre'
        7: cust,  # 'Cust'
        19: race_biwash,  # 'Race Bywash (Cust)'
    }

    for k,v in inflows.items():
        seg_data.loc[k,'flow'] = v
    return seg_data


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


def convert_rch_to_recarray(rch_data): # todo finish when checking
    """
    convert teh reach data to teh record array needed for modflow
    :param rch_data:
    :return:
    """
    raise NotImplementedError


def convert_seg_to_recarray(seg_data): # todo finish when checking
    """
    convert the segment data to the record array needed for modflow
    :return:
    """
    raise NotImplementedError
