# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 14/11/2017 7:34 PM
"""

from __future__ import division
from core import env
from obs_2d_mean_sd import plot_sd_mean_multid, no_change
import socket
import numpy as np
import netCDF4 as nc
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
import pandas as pd
from warnings import warn
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import get_all_well_row_col

def _above_mav(x, **kwargs):
    return x >= 11.3

def _above_half_mav(x, **kwargs):
    return x >= 5.65

def _above_quarter_mav(x, **kwargs):
    return x >= 2.825

def _above_low(x, **kwargs):
    return x >= 1

# todo decide on filters
def plot_all_2d_con(outdir, filter_strs):
    if socket.gethostname() != 'GWATER02':
        raise ValueError('this must be run on GWATER02 as that is where the uncompressed data is stored')

    nc_param_data = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    emma = nc.Dataset(r"C:\mh_waimak_model_data\emma_con.nc")
    bestn = nc.Dataset(r"C:\mh_waimak_model_data\mednload_ucn.nc")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for l in range(smt.layers):
        title = 'coastal component for layer: {:02d}'.format(l)
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=emma,
                                       data_id='coastal', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd': True, 'mean': True, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title))
        plt.close()

    for l in range(smt.layers):
        title = 'inland component for layer: {:02d}'.format(l)
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=emma,
                                       data_id='inland', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd': True, 'mean': True, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title))
        plt.close()

    for l in range(smt.layers):
        title = 'river component for layer: {:02d}'.format(l)
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=emma,
                                       data_id='river', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd': True, 'mean': True, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title))
        plt.close()

    for l in range(smt.layers):
        title = 'best N estimate (g/m3) for layer: {:02d}'.format(l)
        fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                       nc_param_data=nc_param_data, nc_obs_data=bestn,
                                       data_id='inland', function_adjust=no_change,
                                       title=title,
                                       basemap=True, contour={'sd': True, 'mean': True, 'sum': True},
                                       method='mean_sd', contour_color='g')

        fig.savefig(os.path.join(outdir, title))
        plt.close()

    for limit, ftn in zip([1, 2.85, 5.65, 11.3], [_above_low,_above_quarter_mav, _above_half_mav, _above_mav]):
        for l in range(smt.layers):
            title = 'N above {} (g/m3) for layer: {:02d}'.format(limit, l)
            fig, axs = plot_sd_mean_multid(filter_strs=filter_strs, layer=l,
                                           nc_param_data=nc_param_data, nc_obs_data=bestn,
                                           data_id='inland', function_adjust=ftn,
                                           title=title,
                                           basemap=True, contour={'sd': True, 'mean': True, 'sum': True},
                                           method='mean_sd', contour_color='g')

            fig.savefig(os.path.join(outdir, title))
            plt.close()


# both 2D spatial as well as point measurements at wells for example
# take for granted that only 2000 something were run when looking through filters

def plot_well_con(nc_path, outdir, filter_strs): #todo debug
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    plt_data = pd.read_excel(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\WellNSMCPlotNitrateGroups.xlsx")
    all_wells = get_all_well_row_col()
    filter_strs = np.atleast_1d(filter_strs)
    nc_data = nc.Dataset(nc_path)
    supergroups = list(set(plt_data['Supergroup']))

    for supergroup in supergroups:
        fig, axs = plt.subplots(ncols=len(filter_strs), figsize=(18.5, 9.5))
        axs = np.atleast_1d(axs)
        plot_sites = list(plt_data.loc[plt_data['Supergroup'] == supergroup, 'Group'].values)
        well_names = np.array(nc_data.variables['well_name'])
        for filter_str_raw, ax in zip(filter_strs, axs):
            ftype = 0
            filter_str = filter_str_raw
            textadd = ''
            if '~0_' in filter_str_raw:
                filter_str = filter_str_raw.replace('~0_', '')
                ftype = 1
                textadd = 'failed '
            elif '~-1_' in filter_str_raw:
                filter_str = filter_str_raw.replace('~-1', '')
                ftype = 2
                textadd = 'not run '
            elif '~-10_' in filter_str_raw:
                filter_str = filter_str_raw.replace('~-10', '')
                ftype = 3
                textadd = 'not run or failed '

            temp_filter = np.array(nc_data.variables[filter_str])
            if ftype == 0:
                real_filter = temp_filter == 1
            elif ftype == 1:
                real_filter = temp_filter == 0
            elif ftype == 2:
                real_filter = temp_filter == -1
            elif ftype == 3:
                real_filter = temp_filter < 1
            else:
                raise ValueError('shouldnt get here')

            nsmc_nums = np.array(nc_data.variables['nsmc_num'][real_filter])
            all_model_data = []

            # get the data
            labels = []
            for site in plot_sites:
                layer, row, col = all_wells.loc[plt_data.loc[plt_data['Group']==site,'WELL_NO'].values,['layer','row','col']].values.transpose()
                if len(layer) == 0 :
                    warn('no wells for site {} not in observations'.format(site))
                    continue
                labels.append(site)


                model_data = np.array([nc_data.variables['well_obs'][real_filter, l, r, c] for l,r,c in zip(layer, row, col)])
                model_data = model_data.sum(axis=0)
                all_model_data.append(model_data[np.isfinite(model_data)])

            # plot it up
            positions = np.arange(len(labels)) + 1
            t = ax.boxplot(x=all_model_data, positions=positions, labels=labels)
            [[e.set_linewidth(2) for e in j[1]] for j in t.items()]
            ax.set_ylabel('m')
            ax.set_title('{}{}'.format(textadd, filter_str))
        ymax = max([e.get_ylim()[1] for e in axs.flatten()])
        ymin = min([e.get_ylim()[0] for e in axs.flatten()])
        [e.set_ylim(ymin, ymax) for e in axs.flatten()]
        [[tick.set_rotation(45) for tick in ax.get_xticklabels()] for ax in axs.flatten()]

        fig.suptitle(supergroup.title())

        fig.savefig(os.path.join(outdir, supergroup.lower()))
