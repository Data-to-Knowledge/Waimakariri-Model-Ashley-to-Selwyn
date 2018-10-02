# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/09/2018 7:42 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import pandas as pd
import numpy as np


def get_obs_from_nc(data):
    all_vars = data.variables.keys()

    outdata = {}
    outdata['long_name'] = {}
    outdata['target'] = {}
    outdata['std'] = {}
    outdata['sd_type'] = {}
    for var in all_vars:
        if 'well' in var or 'phi' in var:
            continue
        elif data.variables[var].vtype != 'obs':
            continue
        print(var)
        if data.variables[var].units == 'm3/day':
            mult = 1 / 86400
        else:
            mult = 1

        if hasattr(data.variables[var], 'sd_type'):
            outdata['sd_type'][var] = data.variables[var].sd_type
        outdata['long_name'][var] = data.variables[var].long_name
        outdata['std'][var] = 1 / data.variables[var].nsmc_weight * mult  # switch back to sd the weights are 1/sd
        outdata['target'][var] = data.variables[var].target * mult

    outdata = pd.DataFrame(outdata)
    return outdata


def get_param_from_nc(data):
    all_vars = data.variables.keys()

    data_to_collect = ['long_name', 'sd_type', 'initial', 'upper', 'lower', 'p_sd', 'j_sd', ]

    outdata = {}
    for key in data_to_collect:
        outdata[key] = {}
    for var in all_vars:
        if data.variables[var].vtype != 'param':
            continue
        print(var)
        if hasattr(data.variables[var],'units'):
            if data.variables[var].units == 'm3/day':
                mult = 1 / 86400
            else:
                mult = 1

        for dtc in data_to_collect:
            if hasattr(data.variables[var], dtc):
                outdata[dtc][var] = data.variables[var].__getattribute__(dtc)

    outdata = pd.DataFrame(outdata)


    # additional values
    long_vars = ['drn_cond',
                 'rch_mult',
                 'sfr_cond_val',
                 'kv', # decided against
                 'kh',] #decided against

    ids = np.array(data.variables['rch_ppt'])
    for id in ids:
        outdata.loc[id,'j_sd'] = np.nan
    mapper = {key:val for key, val in zip(data.variables['rch_ppt_group'].flag_values,
                                          data.variables['rch_ppt_group'].flag_meanings.split(' '))}
    outdata.loc[ids, 'group'] = np.array([mapper[e] for e in data.variables['rch_ppt_group'][:]])

    group_name = 'rch_ppt'
    outdata.loc[ids, 'initial']= np.array(data.variables['{}_initial'.format(group_name)])
    outdata.loc[ids, 'j_sd']= np.array(data.variables['{}_j_sd'.format(group_name)])
    outdata.loc[ids, 'long_name']= np.array(data.variables['{}'.format(group_name)].long_name)
    outdata.loc[ids, 'lower']= np.array(data.variables['{}_lower'.format(group_name)])
    outdata.loc[ids, 'p_sd']= np.array(data.variables['{}_p_sd'.format(group_name)])
    outdata.loc[ids, 'sd_type']= np.array(data.variables['{}_p_sd'.format(group_name)].sd_type)
    outdata.loc[ids, 'upper']= np.array(data.variables['{}_upper'.format(group_name)])

    group_name = 'drn'
    ids = np.array(data.variables['drns'])
    for id in ids:
        outdata.loc[id,'j_sd'] = np.nan

    outdata.loc[ids, 'initial']= np.array(data.variables['{}_initial'.format(group_name)])
    outdata.loc[ids, 'j_sd']= np.array(data.variables['{}_j_sd'.format(group_name)])
    outdata.loc[ids, 'long_name']= np.array(data.variables['{}_cond'.format(group_name)].long_name)
    outdata.loc[ids, 'lower']= np.array(data.variables['{}_lower'.format(group_name)])
    outdata.loc[ids, 'p_sd']= np.array(data.variables['{}_p_sd'.format(group_name)])
    outdata.loc[ids, 'sd_type']= np.array(data.variables['{}_p_sd'.format(group_name)].sd_type)
    outdata.loc[ids, 'upper']= np.array(data.variables['{}_upper'.format(group_name)])

    group_name = 'sfr'
    ids = np.array(data.variables['sfr_cond'])
    for id in ids:
        outdata.loc[id,'j_sd'] = np.nan

    outdata.loc[ids, 'initial']= np.array(data.variables['{}_initial'.format(group_name)])
    outdata.loc[ids, 'j_sd']= np.array(data.variables['{}_j_sd'.format(group_name)])
    outdata.loc[ids, 'long_name']= np.array(data.variables['{}_cond_val'.format(group_name)].long_name)
    outdata.loc[ids, 'lower']= np.array(data.variables['{}_lower'.format(group_name)])
    outdata.loc[ids, 'p_sd']= np.array(data.variables['{}_p_sd'.format(group_name)])
    outdata.loc[ids, 'sd_type']= np.array(data.variables['{}_p_sd'.format(group_name)].sd_type)
    outdata.loc[ids, 'upper']= np.array(data.variables['{}_upper'.format(group_name)])

    return outdata


if __name__ == '__main__':
    dataset = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    out = get_obs_from_nc(dataset)
    out.to_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\modelling_reports\ashely_waimakarriri_model_build\params_netcdf_obs.csv")
    out2 = get_param_from_nc(dataset)
    out2.to_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\modelling_reports\ashely_waimakarriri_model_build\params_netcdf_params.csv")
