# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 10:25 AM
"""

from __future__ import division
import env
import pandas as pd
import numpy as np
import netCDF4 as nc
import os
import gc
# quick check for the netcdf
def check_missing(nc_path):
    print(nc_path)
    data = nc.Dataset(nc_path)
    variables = list(set(data.variables.keys()) - {u'layer',
                                                   u'row',
                                                   u'col',
                                                   u'nsmc_num',
                                                   u'crs',
                                                   u'latitude',
                                                   u'longitude'})
    nsmc_nums = np.array(data.variables['nsmc_num'])

    outdata = pd.DataFrame(index=nsmc_nums, columns=variables, data=False)
    for var in variables:
        if data[var].ndim == 4:
            temp_data = np.array(data[var][:,5])
        elif data[var].ndim == 3:
            temp_data = np.array(data[var][:])
        else:
            raise ValueError('var: {} has dimensions of {}, why?'.format(var,data[var].ndim))
        for i, num in enumerate(nsmc_nums):
            if i%100 == 0:
                print('{} of {} for {}'.format(i,len(nsmc_nums),var))
            outdata.loc[num,var] = np.isnan(temp_data[i]).all()
        gc.collect()

    outdata.to_csv(os.path.join(os.path.dirname(nc_path),'check_nan_{}.csv'.format(os.path.basename(nc_path).replace('.nc',''))))


if __name__ == '__main__':
    check_missing(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\emma_unc.nc"))
    check_missing(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\post_filter1_hds.nc"))
    check_missing(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\post_filter1_cell_budgets.nc"))