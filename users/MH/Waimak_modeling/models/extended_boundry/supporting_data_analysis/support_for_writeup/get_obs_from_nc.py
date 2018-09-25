# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/09/2018 7:42 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import pandas as pd

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
            mult = 1/86400
        else:
            mult = 1

        if hasattr(data.variables[var],'sd_type'):
            outdata['sd_type'][var] = data.variables[var].sd_type
        outdata['long_name'][var] = data.variables[var].long_name
        outdata['std'][var] = 1/data.variables[var].nsmc_weight*mult  # switch back to sd the weights are 1/sd
        outdata['target'][var] = data.variables[var].target * mult

    outdata = pd.DataFrame(outdata)
    return outdata

if __name__ == '__main__':
    dataset = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    out = get_obs_from_nc(dataset)
    out.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\modelling_reports\ashely_waimakarriri_model_build\params_netcdf_obs.csv")
