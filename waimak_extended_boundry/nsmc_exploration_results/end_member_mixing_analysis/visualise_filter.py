# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 14/11/2017 1:18 PM
"""

from __future__ import division
import env
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
np.random.seed(1)

def plot_river_mixing(outdir, filter_strs):
    filter_strs = np.atleast_1d(filter_strs)
    data = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\nsmc_results\emma_filter\river_end_member_mixing.csv"), index_col=0)

    param_netcdf_file = nc.Dataset(
        env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc"))
    target_data = pd.read_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\4EM_WaimakResults_NSMCFilterWithStreams.csv")
    target_data.loc[:,'Group'] = [e.lower() for e in target_data.Group]
    sites = {}
    for g in set(target_data.SuperGroup):
        sites[g] = list(set(target_data.loc[target_data.SuperGroup==g,'Group']))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for key, ids in sites.items():
        model_data = data.loc[:,ids].transpose()
        fig, axs = plt.subplots(ncols=len(filter_strs), figsize=(18.5, 9.5))
        axs = np.atleast_1d(axs)

        for filter_str_raw, ax in zip(filter_strs, axs.flatten()):
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

            temp_filter = np.array(param_netcdf_file.variables[filter_str])
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
            pos = np.arange(len(ids))+1
            tar_medians = [target_data.loc[target_data.Group==e,'river_med_4EM'].iloc[0] for e in ids]
            tar_mins = [target_data.loc[target_data.Group==e,'river_5th_4EM'].iloc[0] for e in ids]
            tar_maxes = [target_data.loc[target_data.Group==e,'river_95th_4EM'].iloc[0] for e in ids]
            t = ax.boxplot(x=[[l,u] for l,u in zip(tar_mins, tar_maxes)], positions=pos - 0.33, usermedians=tar_medians, whis = 'range',
                           showbox=False, showfliers=False)
            ys = sorted(list(set([e.remove() for e in t['whiskers']])))
            [ax.plot([p-0.33, p-0.33], [l,u], alpha=0.33, color='k') for p, l, u in zip(pos,tar_mins,tar_maxes)]
            [[e.set_alpha(0.33) for e in j[1]] for j in t.items()]

            # plot the models
            usenums = np.in1d(model_data.keys(),np.array(param_netcdf_file.variables['nsmc_num'][real_filter]))
            plt_data = [e[np.isfinite(e)] for e in model_data.loc[:,usenums].values]
            t = ax.boxplot(x=plt_data, positions=pos, labels=[e.title() for e in ids], whis=[5,95])
            [[e.set_linewidth(2) for e in j[1]] for j in t.items()]
            ax.set_ylabel('fraction of alpine river')
            ax.set_title('{}{}'.format(textadd, filter_str))
        fig.suptitle(key.title())
        ymax = max([e.get_ylim()[1] for e in axs.flatten()])
        ymin = min([e.get_ylim()[0] for e in axs.flatten()])
        [[tick.set_rotation(45) for tick in ax.get_xticklabels()] for ax in axs.flatten()]
        [e.set_ylim(ymin, ymax) for e in axs.flatten()]
        fig.savefig(os.path.join(outdir,key))
        plt.close(fig)

if __name__ == '__main__':
    plot_river_mixing(r"T:\Temp\temp_gw_files\test_emma_plots",[u'emma_chch_wt',u'emma_ewf_wt'])