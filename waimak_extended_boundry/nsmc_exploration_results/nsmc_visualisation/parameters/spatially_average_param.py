# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 4:57 PM
"""

from __future__ import division
import env
from waimak_extended_boundry import smt
from waimak_extended_boundry.supporting_data_analysis.NSMC_inputs.recharge_index_array import get_rch_index_array
import numpy as np
from copy import deepcopy
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import netCDF4 as nc
import os

def no_change(x, **kwargs):
    return x

def _log_10(x, **kwargs):
    return np.log10(x)

def extract_data(param_nc, filter_bool, layer, data_id, transform=no_change):
    """
    extract and extrapolate (via cubic spline) data
    :param param_nc: netcdf object for the parameter data
    :param filter_bool: a boolean filter of lenght netchd object
    :param layer: zero indexed layer
    :param data_id: the key for teh netcdf
    :param transform: the function to apply to the interpolated arrays (note that kv/kh is logged prior to this step)
    :return:
    """
    xs, ys = smt.get_model_x_y()
    if data_id in ['kv','kh']:
        kv_x = np.array(param_nc.variables['khv_ppt_x'])
        kv_y = np.array(param_nc.variables['khv_ppt_y'])
        val = np.log10(np.array(param_nc.variables[data_id][filter_bool, layer]))
        temp_mean = np.nanmean(val, axis=0)
        temp_sd = np.nanstd(val, axis=0)
        idx = np.isfinite(temp_mean)
        outmean = griddata(points=(kv_x[idx], kv_y[idx]), values=temp_mean[idx], xi=(xs, ys), method='cubic')
        outsd = griddata(points=(kv_x[idx], kv_y[idx]), values=temp_sd[idx], xi=(xs, ys), method='cubic')

    elif 'rch' in data_id:
        rch_index_array = get_rch_index_array()
        outmean = smt.get_empty_model_grid() * np.nan
        outsd = deepcopy(outmean)
        rchx = np.array(param_nc.variables['rch_ppt_x'])
        rchy = np.array(param_nc.variables['rch_ppt_y'])
        rchgroup = np.array(param_nc.variables['rch_ppt_group'])
        val = np.array(param_nc.variables[data_id][filter_bool])
        for i in set(rchgroup):
            val_idx = np.isclose(i,rchgroup)
            temp_mean = np.nanmean(val[:, val_idx], axis=0)
            temp_sd = np.nanstd(val[:, val_idx], axis=0)
            temp_x = rchx[val_idx]
            temp_y = rchy[val_idx]

            # interpolate with spline and assign data
            grid_idx = np.isclose(rch_index_array,i)
            outmean[grid_idx] = griddata(points=(temp_x, temp_y), values=temp_mean, xi=(xs, ys), method='cubic')[grid_idx]
            outsd[grid_idx] = griddata(points=(temp_x, temp_y), values=temp_sd, xi=(xs, ys), method='cubic')[grid_idx]
    else:
        raise NotImplementedError('{} not yet implemented for 2d plotting'.format(data_id))
    return transform(outmean,layer), transform(outsd,layer)

def plot_sd_mean_multid(filter_strs, layer, nc_param_data, data_id,
                        title, basemap=True, contour={'sd': False, 'mean': False},
                        contour_color='g', vmins=None, vmaxes=None, tranform=_log_10):
    filter_strs = np.atleast_1d(filter_strs)

    # get the data
    plot_data_mean = {}
    plot_data_sd = {}
    textadds = {}
    use_filter_strs = []
    mean_maxes = []
    mean_mins = []
    sd_maxes = []
    sd_mins = []

    for filter_str_raw in filter_strs:
        # id the filter interpretation
        if '~0_' in filter_str_raw:
            filter_str = filter_str_raw.replace('~0_', '')
            ftype = 1
            textadds[filter_str] = 'failed '
        elif '~-1_' in filter_str_raw:
            filter_str = filter_str_raw.replace('~-1', '')
            ftype = 2
            textadds[filter_str] = 'not run '
        elif '~-10_' in filter_str_raw:
            filter_str = filter_str_raw.replace('~-10', '')
            ftype = 3
            textadds[filter_str] = 'not run or failed '
        else:
            ftype = 0
            filter_str = filter_str_raw
            textadds[filter_str] = ''

        # get the filters and handle the possibly different dimensions between the two netcdfs
        temp_filter = np.array(nc_param_data.variables[filter_str])
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
        use_filter_strs.append(filter_str)

        # pull the data out

        mean, sd = extract_data(param_nc=nc_param_data, filter_bool=real_filters, layer=layer, data_id=data_id,transform=tranform)
        plot_data_mean[filter_str] = mean
        plot_data_sd[filter_str] = sd # thought about relative sd, but it makes a mess as the data is zero inflated
        mean_mins.append(np.nanpercentile(mean,1))
        mean_maxes.append(np.nanpercentile(mean,99))
        sd_maxes.append(np.nanpercentile(sd,99))
        sd_mins.append(np.nanpercentile(sd,1))

    # initalised the figure

    ncols = 2
    fig, axs = plt.subplots(ncols=ncols, nrows=len(filter_strs), figsize=(18.5, 9.5))
    axs = np.atleast_2d(axs)
    mean_vmin = min(mean_mins)
    sd_vmin = min(sd_mins)
    if vmins is not None:
        if not isinstance(vmins, dict):
            raise ValueError('vmins must be none or dict')
        try:
            mean_vmin = vmins['mean']
        except KeyError:
            pass
        try:
            sd_vmin = vmins['sd']
        except KeyError:
            pass

    mean_vmax = max(mean_maxes)
    sd_vmax = max(sd_maxes)

    if vmaxes is not None:
        if not isinstance(vmaxes, dict):
            raise ValueError('vmaxes must be none or dict')

        try:
            mean_vmax = vmaxes['mean']
        except KeyError:
            pass
        try:
            sd_vmax = vmaxes['sd']
        except KeyError:
            pass

    for i, fstr in enumerate(use_filter_strs):
        plt_cbar = False
        cbar_lab_sd = None
        cbar_lab_mean = None
        if i == len(use_filter_strs) - 1:
            plt_cbar = True
            cbar_lab_sd = 'Std'
            cbar_lab_mean = 'Mean'

        meanax = axs[i, 0]
        sdax = axs[i, 1]
        smt.plt_matrix(plot_data_mean[fstr], vmin=mean_vmin, vmax=mean_vmax,
                       title='{} mean'.format(fstr),
                       no_flow_layer=layer, ax=meanax, color_bar=plt_cbar, base_map=basemap,
                       cbar_lab=cbar_lab_mean,
                       )

        if contour['mean']:
            x, y = smt.get_model_x_y()
            meanax.contour(x, y, plot_data_mean[fstr], colors=contour_color)

        smt.plt_matrix(plot_data_sd[fstr], vmin=sd_vmin, vmax=sd_vmax,
                       title='{} sd'.format(fstr),
                       no_flow_layer=layer, ax=sdax, color_bar=plt_cbar, base_map=basemap,
                       cbar_lab=cbar_lab_sd,
                       )
        if contour['sd']:
            x, y = smt.get_model_x_y()
            sdax.contour(x, y, plot_data_sd[fstr], colors=contour_color)
        meanax.set_title('Mean {}'.format(fstr))
        sdax.set_title('Std {}'.format(fstr))

    fig.suptitle(title.title())
    fig.tight_layout()
    return fig, axs

def plt_all_spatial_param(outdir,filter_strs):
    nc_param_data = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc"))

    if not os.path.exists(outdir):
        os.makedirs(outdir)


    #kv for each layer
    for l in range(smt.layers):
        title = 'log kv for layer {:02d}'.format(l+1)
        print(title)

        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l, nc_param_data=nc_param_data,
                                       data_id='kv',
                                       title=title, basemap=True, contour={'sd': False, 'mean': False},
                                       contour_color='g', vmins=None, vmaxes=None)
        fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
        plt.close(fig)


        #kh for each layer
    for l in range(smt.layers):
        title = 'log kh for layer {:02d}'.format(l+1)
        print(title)

        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l, nc_param_data=nc_param_data,
                                       data_id='kh',
                                       title=title, basemap=True, contour={'sd': False, 'mean': False},
                                       contour_color='g', vmins=None, vmaxes=None)
        fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
        plt.close(fig)

    # rch
    title = 'recharge multiplier'
    print(title)
    fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=0, nc_param_data=nc_param_data, data_id='rch_mult',
                        title=title, basemap=True, contour={'sd': False, 'mean': False},
                        contour_color='g', vmins=None, vmaxes=None)
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
if __name__ == '__main__':
    plt_all_spatial_param(r"T:\Temp\temp_gw_files\testparam2dplots", filter_strs=['filter2', 'filter3'])