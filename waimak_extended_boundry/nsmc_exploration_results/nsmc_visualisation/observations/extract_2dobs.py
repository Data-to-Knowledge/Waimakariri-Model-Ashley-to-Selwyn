"""
Author: matth
Date Created: 27/11/2017 3:35 PM
"""

from __future__ import division
import numpy as np
import netCDF4 as nc
import os
from waimak_extended_boundry import smt

def extract_2d_obs(filter_str_raw, layer, nc_param_data, nc_obs_data, data_id):
    """

    :param filter_strs: a string or list of strings corresponding to the filter name to use from the param netcdf
                        possible prefixes: ~0_ where the filter failed
                                           ~1_ where the filter was not run
                                           ~10_ where the filter was not run or failed
    :param layer: the zero indexed modflow layer
    :param nc_param_data: the netcdf class of the parameter and observation data
    :param nc_obs_data: the 3-4d netcdf (e.g. concentration or head or drain flux)
    :param data_id: the key for the data
    :param function_adjust: a function to run the data through prior to amalgamation
    :param title: the tile for the plot
    :param basemap: Boolean, argument for ModelTools.plt_matrix
    :param contour: a dictionary keys sum, mean, sd and boolean values of whether or not to contour the data
    :param method: either mean_sd (plots the mean and standard devation as 2 plot)  or sum(plots just eh sum).
    :param contour_color: a string to denote the color for the contouring, default green
    :param vmins: a dictionary with keys sum, mean, sd and keys of the maximum value to use in the pcolor mesh
                  or None(use default min/maxes) not all keys need to present
    :param vmaxes: as vmins, but for the minimum value
    :return: fig,axs (matplotlib figure, array of axes)
    """
    param_nc_nums = np.array(nc_param_data.variables['nsmc_num'])
    obs_nc_nums = np.array(nc_obs_data.variables['nsmc_num'])

    # id the filter interpretation
    if '~0_' in filter_str_raw:
        filter_str = filter_str_raw.replace('~0_', '')
        ftype = 1
    elif '~-1_' in filter_str_raw:
        filter_str = filter_str_raw.replace('~-1', '')
        ftype = 2
    elif '~-10_' in filter_str_raw:
        filter_str = filter_str_raw.replace('~-10', '')
        ftype = 3
    else:
        ftype = 0
        filter_str = filter_str_raw

    # get the filters and handle the possibly different dimensions between the two netcdfs
    temp_filter_full = dict(zip(param_nc_nums, np.array(nc_param_data.variables[filter_str])))
    temp_filter = np.array([temp_filter_full[e] for e in obs_nc_nums])
    if ftype == 0:
        real_filters = temp_filter == 1
    elif ftype == 1:
        real_filters = temp_filter == 0
    elif ftype == 2:
        real_filters = temp_filter == -1
    elif ftype == 3:
        real_filters = temp_filter < 1
    else:
        raise ValueError('shouldnt get here')

    # pull the data out
    if nc_obs_data.variables[data_id].ndim == 4:
        temp = nc_obs_data.variables[data_id][real_filters, layer, :, :]
        if isinstance(temp,np.ma.MaskedArray):
            temp = temp.filled(np.nan)
    elif nc_obs_data.variables[data_id].ndim == 3:
        temp = nc_obs_data.variables[data_id][real_filters]
        if isinstance(temp,np.ma.MaskedArray):
            temp = temp.filled(np.nan)
    else:
        raise ValueError('unexpected data shape')
    return np.nanmean(temp, axis=0)

if __name__ == '__main__':
    nc_obs = nc.Dataset(r"C:\mh_waimak_model_data\mednload_ucn.nc")
    nc_param = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results"
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    for layer in range(smt.layers):
        print(layer+1)
        data = extract_2d_obs('emma_no_wt',layer,nc_param,nc_obs,'mednload')
        smt.array_to_raster(os.path.join(base_dir,'med_n_layer{:02d}.tif'.format(layer+1)),data,layer)