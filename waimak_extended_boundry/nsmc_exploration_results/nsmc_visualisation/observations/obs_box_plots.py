# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 4:59 PM
"""

from __future__ import division
import env
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import os
import pandas as pd
from warnings import warn

np.random.seed(1) # make it so the random distributions don't get mucked up in future

def obs_boxplots(param_netcdf_file, sites, title, yax, filter_strs):
    """

    :param param_netcdf: the parameter netcdf with obs ect (nc object)
    :param sites: an ordered dictionary of lists first list hold everything to plot on the x axis,
                  the second is added together these are string coresponding to the netcdf variable
                  eg. {'eyre': [sfo_e],'eyre2':[sfo_e2],'ohoka':[oho_mill, oho_misc]}
    :param title: title to use
    :param yax: label for the y axis
    :param filter_strs: list of which filter should I use 'str' where multiple filter they appear as column subplots
                        possible prefixes: ~0_ where the filter failed
                                           ~1_ where the filter was not run
                                           ~10_ where the filter was not run or failed
    :return: fig, axs
    """
    filter_strs = np.atleast_1d(filter_strs)
    if not isinstance(sites, OrderedDict):
        raise ValueError('expected ordered dictionary')

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

        nsmc_nums = np.array(param_netcdf_file.variables['nsmc_num'][real_filter])
        all_model_data = []
        all_tar_data = []
        # get the data
        # get the target
        # get the weight/sd
        for site in sites.values():
            model_data = np.zeros(nsmc_nums.shape)
            tar_data = np.zeros(nsmc_nums.shape)
            for subsite in site:
                model_data += np.array(param_netcdf_file[subsite][real_filter])
                u = param_netcdf_file[subsite].target
                wt = param_netcdf_file[subsite].nsmc_weight
                if wt != 0:
                    tar_data += np.random.normal(loc=u, scale=1 / wt, size=len(tar_data))
            if yax == 'm3/s':  # then convert to m3/s from m3/day
                model_data *= 1 / 86400
                tar_data *= 1 / 86400
            all_model_data.append(model_data[np.isfinite(model_data)])
            all_tar_data.append(tar_data)

        # plot the thing
        positions = np.arange(len(all_model_data)) + 1
        # plot the data
        t = ax.boxplot(x=all_tar_data, positions=positions - 0.33, whis=[5,95])
        [[e.set_alpha(0.33) for e in j[1]] for j in t.items()]
        t = ax.boxplot(x=all_model_data, positions=positions, labels=sites.keys(), whis=[5,95])
        [[e.set_linewidth(2) for e in j[1]] for j in t.items()]
        ax.set_ylabel(yax)
        ax.set_title('{}{}'.format(textadd, filter_str))
    ymax = max([e.get_ylim()[1] for e in axs.flatten()])
    ymin = min([e.get_ylim()[0] for e in axs.flatten()])
    [e.set_ylim(ymin, ymax) for e in axs.flatten()]
    fig.suptitle(title.title())

    return fig, axs


def plot_all_obs_boxplots(nc_path, outdir, filter_strs):
    """

    :param nc_path: the path to the param/obs/metadata netcdf
    :param outdir: the directory to save the plots to
    :param filter_strs: a list of filter strings
    :return:
    """
    nc_data = nc.Dataset(nc_path)

    plots = {
        # tile                                  site 1 name         st comps
        'total waimakariri losses': OrderedDict([('Waimakariri', ['sfx_w_all'])]),

        'middle ashley losses': OrderedDict([('Ashley', ['mid_ash_g'])]),

        'silver stream flow at island rd': OrderedDict(
            [('Silver Stream', [u'd_sil_harp', u'd_sil_heyw', u'd_sil_ilnd', ])]),

        'ohoka flow at island rd': OrderedDict([('Ohoka', [u'd_oho_btch', u'd_oho_jefs', u'd_oho_kpoi', u'd_oho_misc',
                                                           u'd_oho_miscz', u'd_oho_mlbk', u'd_oho_whit', ])]),

        'cust flow at threlkelds rd': OrderedDict([('Cust', [u'sfo_c_tlks'])]),

        'cam at youngs rd': OrderedDict([('Cam', [u'd_cam_mrsh', u'd_cam_yng', u'd_nbk_mrsh', u'd_sbk_mrsh'])]),
        'taranaki at stokes': OrderedDict([('Taranaki', [u'd_tar_gre', u'd_tar_stok', ])]),

        'waikuku at leggitts': OrderedDict([('Waikuku', [u'd_kuku_leg', ])]),

        'saltwater at factory rd': OrderedDict([('Saltwater', [u'd_salt_fct', u'd_salt_s', u'd_salt_top'])]),

        'total selwyn gains': OrderedDict([('Selwyn Steams', [u'sel_str'])]),

        'eyre river targets': OrderedDict([('Wollffs Rd', [u'sfo_e_wolf']),
                                           ('Poyntzs Rd', ['sfo_e_poyn']),
                                           ('Downs Rd', ['sfo_e_down']),
                                           ('S. Eyre Rd', ['sfo_e_seyr'])]),

        'Waimakariri Reaches': OrderedDict([('Courtenay', [u'sfx_w1_cou']),
                                            ('Halkett', ['sfx_w2_tom']),
                                            ('Weedons Ross Rd', ['sfx_w3_ros']),
                                            ('Crossbank', ['sfx_w4_mcl']),
                                            ('Wrights Cut', ['sfx_w5_wat']),
                                            ('Old SH1', ['sfx_w6_sh1'])]),

        'Carpet Drains': OrderedDict([('Ashley Carpet', ['d_ash_c']),
                                      ('Cust Carpet', ['d_cust_c'])]),

        'Total Chch Flow': OrderedDict([('Chch Streams', ['chch_str'])]),

        'Coastal Discharge': OrderedDict([('Ashley', ['chb_ash']),
                                          ('Cust', ['chb_cust']),
                                          ('Chch', ['chb_chch']),
                                          ('Selwyn', ['sel_off'])]),

        'Inland Vertical Targets': OrderedDict([(e, [e]) for e in [u'brnthl_4_1',
                                                                   u'eyrftn_2_1',
                                                                   u'eyrftn_6_2',
                                                                   u'eyrftn_6_1',
                                                                   u'eyrftm_2_1',
                                                                   u'eyrftl_7_2',
                                                                   u'oxfds_6_2',
                                                                   u'oxfdnr_3_1',
                                                                   u'oxfdnr_4_3',
                                                                   u'oxfdnr_4_1',
                                                                   ]]),

        'Coastal Vertical Targets': OrderedDict([(e, [e]) for e in [u'wdnd_4_2',
                                                                    u'wdnd_8_4',
                                                                    u'wdnd_8_2',
                                                                    u'peg_8_7',
                                                                    u'peg_9_8',
                                                                    u'peg_10_9',
                                                                    u'peg_10_7',
                                                                    u'chch_4_2'
                                                                    ]]),

        'Waimakariri Springfed Streams': OrderedDict([('Saltwater', [u'd_salt_fct', u'd_salt_s', u'd_salt_top']),
                                                      ('Waikuku', [u'd_kuku_leg', ]),
                                                      ('Taranaki', [u'd_tar_gre', u'd_tar_stok', ]),
                                                      ('Cam',
                                                       [u'd_cam_mrsh', u'd_cam_yng', u'd_nbk_mrsh', 'd_sbk_mrsh']),
                                                      ('Ohoka',
                                                       [u'd_oho_btch', u'd_oho_jefs', u'd_oho_kpoi', u'd_oho_misc',
                                                        u'd_oho_miscz', u'd_oho_mlbk', u'd_oho_whit', ]),
                                                      ('Silver Stream', [u'd_sil_harp', u'd_sil_heyw', u'd_sil_ilnd', ])
                                                      ]),
    }

    for key, val in plots.items():
        if key in ['Coastal Vertical Targets', 'Inland Vertical Targets']:
            ylab = 'm'
        else:
            ylab = 'm3/s'

        fig, axs = obs_boxplots(param_netcdf_file=nc_data,
                                sites=val, title=key, yax=ylab,
                                filter_strs=filter_strs)

        # save figures
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        fig.savefig(os.path.join(outdir, fig._suptitle._text+'.png'))
        plt.close(fig)


def plot_hds_boxplots(nc_path, outdir, filter_strs):
    """
    plot up heads boxplots
    :param nc_path: the param/obs ect netcdf
    :param outdir: directory to plot
    :param filter_strs: strings to use
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    filter_strs = np.atleast_1d(filter_strs)
    nc_data = nc.Dataset(nc_path)
    plt_data = pd.read_excel(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\WellNSMCPlotGroups.xlsx")
    groups = list(set(plt_data['Group_']))

    for group in groups:
        fig, axs = plt.subplots(ncols=len(filter_strs), figsize=(18.5, 9.5))
        axs = np.atleast_1d(axs)
        plot_sites = list(plt_data.loc[plt_data['Group_'] == group, 'Well'].values)
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
            all_target_data = []

            # get the data
            labels = []
            for subsite in plot_sites:
                well_idx = np.where(well_names==subsite)[0]
                if len(well_idx) == 0:
                    warn('well {} not in observations'.format(subsite))
                    continue
                if len(well_idx) != 1:
                    raise ValueError('multiple well indexing is not supported for {}'.format(subsite))

                labels.append(subsite)
                well_idx = well_idx[0]

                tar_data = np.zeros(nsmc_nums.shape)
                model_data = np.zeros(nsmc_nums.shape)
                model_data += np.array(nc_data.variables['well_obs'][real_filter,well_idx])
                u = nc_data['well_target'][well_idx]
                wt = nc_data['well_weight'][well_idx]
                if wt != 0:
                    tar_data += np.random.normal(loc=u, scale=1 / wt, size=len(tar_data))
                all_model_data.append(model_data[np.isfinite(model_data)])
                all_target_data.append(tar_data)

            # plot it up
            positions = np.arange(len(labels)) + 1
            t = ax.boxplot(x=all_target_data, positions=positions - 0.33, whis=[5,95])
            [[e.set_alpha(0.33) for e in j[1]] for j in t.items()]
            t = ax.boxplot(x=all_model_data, positions=positions, labels=labels, whis=[5,95])
            [[e.set_linewidth(2) for e in j[1]] for j in t.items()]
            ax.set_ylabel('m')
            ax.set_title('{}{}'.format(textadd, filter_str))
        ymax = max([e.get_ylim()[1] for e in axs.flatten()])
        ymin = min([e.get_ylim()[0] for e in axs.flatten()])
        [e.set_ylim(ymin, ymax) for e in axs.flatten()]
        [[tick.set_rotation(45) for tick in ax.get_xticklabels()] for ax in axs.flatten()]

        fig.suptitle(group.title())

        fig.savefig(os.path.join(outdir, group.lower()+'.png'))
        plt.close(fig)



def plot_obs_all_all_boxplots(outdir, filters):
    """
    quick wrapper to quicly let me plot everything
    :param outdir:
    :param filters:
    :return:
    """
    nc_path = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    plot_all_obs_boxplots(nc_path,outdir,filters)
    plot_hds_boxplots(nc_path,outdir,filters)


if __name__ == '__main__':
    plot_obs_all_all_boxplots(r"T:\Temp\temp_gw_files\testalllobsbox",['emma_ewf_wt', 'run_mt3d'])