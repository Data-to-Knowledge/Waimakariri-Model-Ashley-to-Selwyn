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


def get_base_sfr_data():
    """
    get the pre-optimisation data, exlcude hcond etc????
    :return:
    """
    seg_data = pd.read_hdf(os.path.join(env.sdp_required, 'sfr_data.hdf'), 'seg_data')
    rch_data = pd.read_hdf(os.path.join(env.sdp_required, 'sfr_data.hdf'), 'reach_data')

    return seg_data, rch_data


def add_sfr_cond(model_id, rch_data):
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
    cond_dict = {k: v for k, v in zip(conduct_names, conduct_vals)}

    # get the parameter mapper
    param_mapper = pd.read_table(os.path.join(env.sdp_required, 'base_for_nsmc_real/sfr_segdata.tpl'),
                                 skiprows=1, index_col=0)

    for nseg, hcond1_nm, hcond2_nm in param_mapper[['hcond1', 'hcond2']].itertuples(True, None):
        hcond1 = cond_dict[hcond1_nm.replace('$', '').lower()]
        hcond2 = cond_dict[hcond2_nm.replace('$', '').lower()]
        m = (hcond2 - hcond1) / rch_data.loc[rch_data.iseg == nseg, 'rchlen'].sum()  # slope of conductance down segment
        # Series of conductance linearly varying from upper to lower
        rhcond = (m * (rch_data.loc[rch_data.iseg == nseg, 'rchlen'].cumsum() -
                       rch_data.loc[rch_data.iseg == nseg, 'rchlen'] / 2.0)) \
                 + hcond1  # m*x+c - x is distance along segment to cell midpoit

        rch_data.loc[rch_data.iseg == nseg, 'strhc1'] = rhcond

    return rch_data


def add_sfr_inflows(model_id, seg_data):
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
    race_biwash = float(param_data.variables['mid_c_flo'][nidx])

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

    for k, v in inflows.items():
        seg_data.loc[k, 'flow'] = v
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


def convert_rch_to_recarray(rch_data):
    """
    convert teh reach data to teh record array needed for modflow
    :param rch_data:
    :return:
    """
    dtype = [('node', '<i4'), ('k', '<i4'), ('i', '<i4'), ('j', '<i4'), ('iseg', '<i4'), ('ireach', '<i4'),
             ('rchlen', '<f4'), ('strtop', '<f4'), ('slope', '<f4'), ('strthick', '<f4'), ('strhc1', '<f4'),
             ('thts', '<f4'), ('thti', '<f4'), ('eps', '<f4'), ('uhc', '<f4'), ('reachID', '<i4'), ('outreach', '<i4')]
    keys = [k[0] for k in dtype]
    dtype = {k[0]: k[1] for k in dtype}
    rch_data = rch_data.loc[:, keys].to_records(index=False, column_dtypes=dtype)
    return rch_data


def convert_seg_to_recarray(seg_data):
    """
    convert the segment data to the record array needed for modflow
    :return:
    """
    dtype = [('nseg', '<i4'), ('icalc', '<i4'), ('outseg', '<i4'), ('iupseg', '<i4'), ('iprior', '<i4'),
             ('nstrpts', '<i4'),
             ('flow', '<f4'), ('runoff', '<f4'), ('etsw', '<f4'), ('pptsw', '<f4'), ('roughch', '<f4'),
             ('roughbk', '<f4'),
             ('cdpth', '<f4'), ('fdpth', '<f4'), ('awdth', '<f4'), ('bwdth', '<f4'), ('hcond1', '<f4'),
             ('thickm1', '<f4'),
             ('elevup', '<f4'), ('width1', '<f4'), ('depth1', '<f4'), ('thts1', '<f4'), ('thti1', '<f4'),
             ('eps1', '<f4'),
             ('uhc1', '<f4'), ('hcond2', '<f4'), ('thickm2', '<f4'), ('elevdn', '<f4'), ('width2', '<f4'),
             ('depth2', '<f4'),
             ('thts2', '<f4'), ('thti2', '<f4'), ('eps2', '<f4'), ('uhc2', '<f4')]
    keys = [k[0] for k in dtype]
    dtype = {k[0]: k[1] for k in dtype}
    seg_data.reset_index(inplace=True)
    seg_data = seg_data.loc[:, keys].to_records(index=False, column_dtypes=dtype)

    return seg_data


if __name__ == '__main__':

    print 'done'
