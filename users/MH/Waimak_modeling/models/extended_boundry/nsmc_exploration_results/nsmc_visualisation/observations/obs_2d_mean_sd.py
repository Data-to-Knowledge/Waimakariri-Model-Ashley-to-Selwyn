# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 5:07 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import os
import socket

def no_change(x, **kwargs):
    """
    small fuction to add to transformations
    :param x:
    :param kwargs: allows other kwargs to showup
    :return: x unchnaged
    """
    return x


def is_dry(x, layer):
    """
    small function, Returns boolean array True if cell is dry
    :param x:
    :param layer:
    :return:
    """
    no_flow = smt.get_no_flow(layer)
    no_flow[no_flow < 0] = 0
    elv_db = smt.calc_elv_db()
    temp = (x<elv_db[layer+1]) & no_flow.astype(bool)
    return temp


def neg_zero_pos(x, **kwargs):
    """
    small function returns array of -1, 0, 1 for negative zero and positive values respectivly
    :param x: array
    :param kwargs: allows things to be passed to other functions
    :return:
    """
    lower_lim = 0
    upper_lim = 0
    outdata = np.zeros(x.shape, np.int8)
    outdata[x < lower_lim] = -1
    outdata[x > upper_lim] = 1
    return outdata.astype(np.float16)


def plot_sd_mean_multid(filter_strs, layer, nc_param_data, nc_obs_data, data_id, function_adjust,
                        title, basemap=True, contour={'sd': True, 'mean': True, 'sum': True}, method='mean_sd',
                        contour_color='g', vmins=None, vmaxes=None):
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
    filter_strs = np.atleast_1d(filter_strs)
    param_nc_nums = np.array(nc_param_data.variables['nsmc_num'])
    obs_nc_nums = np.array(nc_obs_data.variables['nsmc_num'])

    # get the data
    plot_data_mean = {}
    plot_data_sd = {}
    plot_data_sum = {}
    textadds = {}
    use_filter_strs = []
    mean_maxes = []
    mean_mins = []
    sd_maxes = []
    sd_mins = []
    sum_maxes = []
    sum_mins = []

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
        use_filter_strs.append(filter_str)

        # pull the data out
        if nc_obs_data.variables[data_id].ndim == 4:
            temp = nc_obs_data.variables[data_id][real_filters, layer, :, :]
            if isinstance(temp,np.ma.MaskedArray):
                temp = temp.filled(np.nan)
            temp = function_adjust(temp, layer=layer)
        elif nc_obs_data.variables[data_id].ndim == 3:
            temp = nc_obs_data.variables[data_id][real_filters]
            if isinstance(temp,np.ma.MaskedArray):
                temp = temp.filled(np.nan)
            temp = function_adjust(temp, layer=0)
        else:
            raise ValueError('{} does not have 3 or 4 dims'.format(data_id))

        mean = np.nanmean(temp, axis=0)
        sd = np.nanstd(temp, axis=0) # thought about relative sd, but the data is zero inflated so the scale gets buggered
        dsum = np.nansum(temp, axis=0)
        plot_data_mean[filter_str] = mean
        plot_data_sd[filter_str] = sd
        plot_data_sum[filter_str] = dsum
        sum_maxes.append(np.nanpercentile(dsum,99))
        sum_mins.append(np.nanpercentile(dsum,1))
        mean_maxes.append(np.nanpercentile(mean,99))
        mean_mins.append(np.nanpercentile(mean,1))
        sd_maxes.append(np.nanpercentile(sd,99))
        sd_mins.append(np.nanpercentile(sd,1))

    # initalised teh figure
    if method == 'mean_sd':
        ncols = 2
    elif method == 'sum':
        ncols = 1
    else:
        raise ValueError('unexpected value for method')
    # initalise the plotting options
    fig, axs = plt.subplots(ncols=ncols, nrows=len(filter_strs), figsize=(18.5, 9.5))
    axs = np.atleast_2d(axs)
    mean_vmin = min(mean_mins)
    sd_vmin = min(sd_mins)
    sum_vmin = min(sum_mins)

    # get any user specified vmin/maxes
    if vmins is not None:
        if not isinstance(vmins,dict):
            raise ValueError('vmins must be none or dict')
        try:
           mean_vmin = vmins['mean']
        except KeyError:
            pass
        try:
           sd_vmin = vmins['sd']
        except KeyError:
            pass
        try:
           sum_vmin = vmins['sum']
        except KeyError:
            pass

    mean_vmax = max(mean_maxes)
    sd_vmax = max(sd_maxes)
    sum_vmax = max(sum_maxes)

    if vmaxes is not None:
        if not isinstance(vmaxes,dict):
            raise ValueError('vmaxes must be none or dict')

        try:
           mean_vmax = vmaxes['mean']
        except KeyError:
            pass
        try:
           sd_vmax = vmaxes['sd']
        except KeyError:
            pass
        try:
           sum_vmax = vmaxes['sum']
        except KeyError:
            pass

    for i, fstr in enumerate(use_filter_strs):
        # no color bar on the all but the last set
        plt_cbar = False
        cbar_lab_sd = None
        cbar_lab_mean = None
        cbar_lab_sum = None
        if i == len(use_filter_strs) - 1:  # set up the color bar at the end of each set
            plt_cbar = True
            cbar_lab_sd = 'Std'
            cbar_lab_mean = 'Mean'
            cbar_lab_sum = 'Sum'

        # plot the data
        if method == 'mean_sd':
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
        elif method == 'sum':
            sumax = axs[0, i]
            smt.plt_matrix(plot_data_sum[fstr], vmin=sum_vmin, vmax=sum_vmax,
                           title='{} sum'.format(fstr),
                           no_flow_layer=layer, ax=sumax, color_bar=plt_cbar, base_map=basemap,
                           cbar_lab=cbar_lab_sum,
                           )

            if contour['sum']:
                x, y = smt.get_model_x_y()
                sumax.contour(x, y, plot_data_sum[fstr], colors=contour_color)
            sumax.set_title('sum {}'.format(fstr))
        else:
            raise ValueError("shouldn't get here")

    fig.suptitle(title.title())
    fig.tight_layout()
    return fig, axs


def plot_all_2d_obs(outdir,filter_strs):
    """
    plot all of the 2d observations for this set of nsmc
    :param outdir: the directory to save the plots in (name is title.png)
    :param filter_strs: list of filter strings
                        possible prefixes: ~0_ where the filter failed
                                           ~1_ where the filter was not run
                                           ~10_ where the filter was not run or failed
    :return:
    """
    if socket.gethostname() != 'GWATER02':
        raise ValueError('this must be run on GWATER02 as that is where the uncompressed data is stored')
        # compressed data takes way to long to query

    nc_param_data = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    nc_hds_data = nc.Dataset(r"C:\mh_waimak_model_data\post_filter1_hds.nc")
    nc_bud_data = nc.Dataset(r"C:\mh_waimak_model_data\post_filter1_budget.nc")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # hds by layer
    for l in range(smt.layers):
        title = 'heads for layer {:02d}'.format(l+1)
        print('plotting {}'.format(title))
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=nc_hds_data,
                                       data_id='heads', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd':False, 'mean': True, 'sum': True},
                                       method='mean_sd', contour_color='g')


        fig.savefig(os.path.join(outdir, title.replace(' ','_')+'.png'))
        plt.close(fig)

    # drn flux
    if False: # quickly stop from running
        title = 'drain flux'
        print('plotting {}'.format(title))
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=0,
                                       nc_param_data=nc_param_data, nc_obs_data=nc_bud_data,
                                       data_id='drains', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd': False, 'mean': False, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title.replace(' ','_')+'.png'))
        plt.close(fig)

    # layer one bottom flux (as 1 0 -1)
    for l in range(smt.layers):
        title = 'layer {:02d} up vs downflow'.format(l+1)
        print('plotting {}'.format(title))
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=nc_bud_data,
                                       data_id='FLOW LOWER FACE'.lower(), function_adjust=neg_zero_pos,
                                       title=title,
                                       basemap=True, contour={'sd': False, 'mean': False, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title.replace(' ','_')+'.png'))
        plt.close(fig)

    # where the model is dry
    for l in range(smt.layers):
        title = 'drycells for layer {:02d}'.format(l+1)
        print('plotting {}'.format(title))
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=nc_hds_data,
                                       data_id='heads', function_adjust=is_dry,
                                       title=title,
                                       basemap=True, contour={'sd': False, 'mean': False, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title.replace(' ','_')+'.png'))
        plt.close(fig)


if __name__ == '__main__':
    plot_all_2d_obs(r"T:\Temp\temp_gw_files\test2dplots",    filter_strs = ['emma_converge', '~0_emma_converge'])