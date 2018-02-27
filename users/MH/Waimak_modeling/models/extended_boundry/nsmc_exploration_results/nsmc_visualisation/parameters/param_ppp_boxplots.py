# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 3:28 PM
"""

from __future__ import division
from core import env
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

np.random.seed(1) # stop the distributions from changing between iterations

def _no_change(x):  # note that the transformations act on lists of numpy arrays
    return x

def _to_m3s(x):
    """
    convert to m3/s from m3/day
    :param x:
    :return:
    """
    x = np.atleast_2d(x)
    return [e for e in x/86400]
def _log10(x):
    """
    log transform x
    :param x:
    :return:
    """
    x = np.atleast_2d(x)
    return [e for e in np.log10(x)]


def gen_dist(sd_type, u, sd, n):
    """
    generate teh appropriate distribution
    :param sd_type: 'lin' or 'log'
    :param u: mean (not log transformed)
    :param sd: standard deviation (possibly log transformed
    :param n: number of samplse to create
    :return: randomlly distributed parameters
    """
    if sd_type == 'log':
        out = np.random.lognormal(np.log(u), sd/np.log10(np.e), size=n)
    elif sd_type == 'lin':
        out = np.random.normal(u, sd, size=n)
    else:
        raise ValueError('{} is unexpected sd_type'.format(sd_type))
    return out


def _plt_parm_boxplot(ax, opt_prior, post_opt_dist, post_jac_dist, labels, filter_dist, limits, ylab, ax_title):
    """
    plot up the parameter boxplots
    :param ax: the matplotlib ax object
    :param opt_prior: the optimisation prior
    :param post_opt_dist: list of numpy arrays the distribution for the priors
    :param post_jac_dist: list of numpy arrays the distibution after the sensitivitiy analysis
    :param labels: the labels to plot list
    :param filter_dist: the data for the distribution after the filtering
    :param limits: list of the minium and maximum [[min1,max1], [min2,max2],...[minn,maxn]]
    :param ylab: the string to put on the ylabel
    :param ax_title: the title for the axis
    :return:
    """
    # watch inputs carefully
    if not all([len(e) == len(opt_prior) for e in [opt_prior, post_opt_dist, post_jac_dist,
                                                   filter_dist, labels, limits]]):
        raise ValueError('all distributions must have the same length')

    # the simple plotting set up
    stp=1.5
    positions = np.arange(0,len(opt_prior)*stp, stp) + 1
    # plot the data
    # post_opt_dist
    t = ax.boxplot(x=post_opt_dist, positions=positions - 0.33, whis=[5,95])
    [[e.set_alpha(0.33) for e in j[1]] for j in t.items()]

    # post_jac_dist
    t = ax.boxplot(x=post_jac_dist, positions=positions, whis=[5,95])
    [[e.set_alpha(0.33) for e in j[1]] for j in t.items()]

    # filter_dist
    t = ax.boxplot(x=filter_dist, positions=positions + 0.33, whis=[5,95])
    [[e.set_linewidth(2) for e in j[1]] for j in t.items()]

    # plot limits
    for lim, pos in zip(limits, positions):
        # plot upper
        ax.plot([pos - 0.33, pos + 0.33], [max(np.atleast_1d(lim)), max(np.atleast_1d(lim))], color='k', linewidth=2)
        # plot lower
        ax.plot([pos - 0.33, pos + 0.33], [min(np.atleast_1d(lim)), min(np.atleast_1d(lim))], color='k', linewidth=2)

    # opt_prior
    t = ax.boxplot(x=opt_prior, positions=positions, usermedians=opt_prior, labels=labels,
                   showcaps=False, showbox=False, showfliers=False, medianprops={'color': 'm', 'linewidth': 2})

    ax.set_xlim([0, positions[-1] + 1])
    ymax = np.max(limits)
    ymin = np.min(limits)
    ax.set_ylim(ymin - 0.05 * np.abs(ymin), ymax + 0.05 * np.abs(ymax))


    ax.set_ylabel(ylab)
    ax.set_title(ax_title)


def _get_plting_dists(nc_file, sites, data_types, filter_array, layers):
    """

    :param nc_file: param netcdf
    :param sites: dictionary {site: [components]}
    :param data_types: array of one of rch, kv, kh, sfr, drn, simple
    :return: opt_prior, post_opt_dist, post_jac_dist, filter_dist, labels, limits
    """
    opt_prior, post_opt_dist, post_jac_dist, filter_dist, labels, limits = [], [], [], [], [], []
    sort_idx = np.argsort(sites.keys())
    for data_type, site, layer in zip(np.array(data_types)[sort_idx], np.array(sites.keys())[sort_idx],
                                      np.array(layers)[sort_idx]):
        labels.append(site)
        array_len = filter_array.sum()
        opt_prior_t = np.array([0])[np.newaxis, :]
        limits_t = np.array([0, 0])[np.newaxis, :]
        post_opt_dist_t = np.zeros((array_len))[np.newaxis, :]
        post_jac_dist_t = np.zeros((array_len))[np.newaxis, :]
        filter_dist_t = np.zeros((array_len))[np.newaxis, :]
        for feature in sites[site]:
            # for each feature I need
            if data_type == 'simple':
                sd_type = nc_file[feature].sd_type
                opt_p = nc_file[feature].opt_p
                p_sd = nc_file[feature].p_sd
                j_sd = nc_file[feature].j_sd
                u = nc_file[feature].initial
                up = nc_file[feature].upper
                lo = nc_file[feature].lower
                addju = 0
                if feature == 'llrzf':
                    u +=1000000
                    addju = - 1000000

                opt_prior_t = np.concatenate((opt_prior_t, [[opt_p]]), axis=0)  # axis just to show for future
                limits_t = np.concatenate((limits_t, [[lo, up]]))
                post_opt_dist_t = np.concatenate(
                    (post_opt_dist_t, gen_dist(sd_type, u, p_sd, array_len)[np.newaxis, :]+addju))
                post_jac_dist_t = np.concatenate(
                    (post_jac_dist_t, gen_dist(sd_type, u, j_sd, array_len)[np.newaxis, :]+addju))
                filter_dist_t = np.concatenate((filter_dist_t, np.array(nc_file[feature][filter_array])[np.newaxis,:]))

            else:
                if data_type == 'rch':
                    id_str = 'rch_ppt'
                    prefix = 'rch_ppt'
                    val_str = 'rch_mult'

                elif data_type == 'kv':
                    id_str = 'khv_ppt'
                    prefix = val_str = 'kv'

                elif data_type == 'kh':
                    id_str = 'khv_ppt'
                    prefix = val_str = 'kh'

                elif data_type == 'sfr':
                    # sfr_lower
                    id_str = 'sfr_cond'
                    prefix = 'sfr'
                    val_str = 'sfr_cond_val'

                elif data_type == 'drn':  # I'm not going to implement now
                    id_str = 'drns'
                    prefix = 'drn'
                    val_str = 'drn_cond'
                else:
                    raise ValueError('wrong input for data_type {}'.format(data_type))
                if data_type in ['kv', 'kh']:
                    ids = np.array(nc_file[id_str])
                    idx = ids == feature
                    sd_type = nc_file['{}_p_sd'.format(prefix)].sd_type
                    opt_p = nc_file['{}_opt_p'.format(prefix)][layer][idx][0]
                    p_sd = nc_file['{}_p_sd'.format(prefix)][layer][idx][0]
                    j_sd = nc_file['{}_j_sd'.format(prefix)][layer][idx][0]
                    u = nc_file['{}_initial'.format(prefix)][layer][idx][0]
                    up = nc_file['{}_upper'.format(prefix)][layer][idx][0]
                    lo = nc_file['{}_lower'.format(prefix)][layer][idx][0]
                    filter_dist_t = np.concatenate(
                        (filter_dist_t, np.array(nc_file[prefix][filter_array, layer, idx]).transpose()))
                else:
                    ids = np.array(nc_file[id_str])
                    idx = ids == feature
                    sd_type = nc_file['{}_p_sd'.format(prefix)].sd_type
                    opt_p = nc_file['{}_opt_p'.format(prefix)][idx][0]
                    p_sd = nc_file['{}_p_sd'.format(prefix)][idx][0]
                    j_sd = nc_file['{}_j_sd'.format(prefix)][idx][0]
                    u = nc_file['{}_initial'.format(prefix)][idx][0]
                    up = nc_file['{}_upper'.format(prefix)][idx][0]
                    lo = nc_file['{}_lower'.format(prefix)][idx][0]
                    filter_dist_t = np.concatenate((filter_dist_t, np.array(nc_file['rch_mult'][filter_array, idx]).transpose()))

                opt_prior_t = np.concatenate((opt_prior_t, [[opt_p]]), axis=0)  # axis just to show for future
                limits_t = np.concatenate((limits_t, [[lo, up]]))
                post_opt_dist_t = np.concatenate(
                    (post_opt_dist_t, gen_dist(sd_type, u, p_sd, array_len)[np.newaxis, :]))
                post_jac_dist_t = np.concatenate(
                    (post_jac_dist_t, gen_dist(sd_type, u, j_sd, array_len)[np.newaxis, :]))

        opt_prior.append(np.nanmean(opt_prior_t[1:], axis=0))
        limits.append(np.nanmean(limits_t[1:], axis=0))
        post_opt_dist.append(np.nanmean(post_opt_dist_t[1:], axis=0))
        post_jac_dist.append(np.nanmean(post_jac_dist_t[1:], axis=0))
        filter_dist.append(np.nanmean(filter_dist_t[1:], axis=0))

    return opt_prior, post_opt_dist, post_jac_dist, filter_dist, labels, limits


def plot_supergroup_ppp_boxplots(filter_strs, param_nc, layers, supergroup, sites, datatypes, ylab,
                                 transformation=_no_change):
    """

    :param filter_strs: typical filter string arguments
    :param param_nc: the netcdf object for the parameter file
    :param layers: list of len sites.keys{} layer to use (assumed to be constant)
    :param supergroup: name of the super group used for the title
    :param sites: the sites {name: [site names to amagamate]}
    :param datatypes: list of len sites.keys() (assumed ot be one value
    :param ylab: label to put on data
    :param transformation: the transformation to do to the distributions (E.g. log10) (a function)
    :return: fig, axs
    """
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

        temp_filter = np.array(param_nc.variables[filter_str])
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

        # get the data
        opt_prior, post_opt_dist, post_jac_dist, \
        filter_dist, labels, limits = _get_plting_dists(nc_file=param_nc,
                                                        sites=sites,
                                                        data_types=datatypes,
                                                        filter_array=real_filter,
                                                        layers=layers)

        # apply transformation
        opt_prior = transformation(opt_prior)
        post_opt_dist= transformation(post_opt_dist)
        post_jac_dist= transformation(post_jac_dist)
        filter_dist= transformation(filter_dist)
        limits= transformation(limits)

        _plt_parm_boxplot(ax, opt_prior, post_opt_dist, post_jac_dist, labels, filter_dist, limits, ylab,
                          ax_title='{}{}'.format(textadd, filter_str))

        ymax = max([e.get_ylim()[1] for e in axs.flatten()])
        ymin = min([e.get_ylim()[0] for e in axs.flatten()])
        [e.set_ylim(ymin, ymax)
         for e in axs.flatten()]
    fig.suptitle(supergroup.title())
    [[tick.set_rotation(45) for tick in ax.get_xticklabels()] for ax in axs.flatten()]

    return fig, axs


def plot_all_ppp_boxplots(outdir, filter_strs):
    """
    plot up all of the ppp boxplots
    :param outdir: dir to save the files in
    :param filter_strs: typical filter strings argument
    :return:
    """
    filter_strs = np.atleast_1d(filter_strs)
    khv_rch_points = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\LSR_KhKv_PilotPoints_Groups.xlsx')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    param_nc = nc.Dataset(env.gw_met_data("mh_modeling/netcdfs_of_key_modeling_data/nsmc_params_obs_metadata.nc"))

    # southern boundry fluxes
    title = 'SW Boundary Flux'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0,0], supergroup=title,
                                           sites={'Inland':['ulrzf'], 'Coastal':['llrzf']},
                                           datatypes = ['simple', 'simple'], ylab='m3/s',
                                           transformation=_to_m3s)
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)

    # PUMPING
    title = 'pumping multipliers'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0,0,0], supergroup=title,
                                           sites={'CHCH':['pump_c'], 'Waimakariri':['pump_w'], 'Selwyn':['pump_s']},
                                           datatypes = ['simple', 'simple', 'simple'], ylab='percent pumping')
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
    # nraces , nboundry flux, s river, south races
    title = 'other multipliers'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0,0,0,0], supergroup=title,
                                           sites={'WIL Races':['n_race'], 'N Boundary':['nbndf'],
                                                  'Selwyn Streams':['sriv'], 'Selwyn Races':['s_race']},
                                           datatypes = ['simple', 'simple', 'simple', 'simple'],
                                           ylab='percent injection')
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)

    # river inflows (may need to split due to scale)
    title = 'Eyre Flow'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0], supergroup=title,
                                           sites={'Eyre Inflow':['top_e_flo']},
                                           datatypes = ['simple'], ylab='m3/s',
                                           transformation=_to_m3s)
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
    title = 'Cust Inflows'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0,0], supergroup=title,
                                           sites={'Top Cust':['top_c_flo'], 'WIL Biwash':['mid_c_flo']},
                                           datatypes = ['simple', 'simple'], ylab='m3/s',
                                           transformation=_to_m3s)
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
    # fault multipliers
    title = 'Log Fault Multipliers'
    print(title)
    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers = [0,0], supergroup=title,
                                           sites={'Kh':['fkh_mult'], 'Kv':['fkv_mult']},
                                           datatypes = ['simple', 'simple'], ylab='log of multiplier',
                                           transformation=_log10)
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
    # recharge (tricky)
    title = 'Land Surface Recharge Multiplier'
    print(title)
    sites = {}
    for g in set(khv_rch_points.loc[khv_rch_points['dtype']=='rch','Group']):
        sites[g] = list(set(khv_rch_points.loc[khv_rch_points['Group']==g,'id']))
    layers = np.zeros(len(sites.keys()))
    datatypes = ['rch' for e in sites.keys()]

    fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                           layers=layers, supergroup=title,
                                           sites=sites,
                                           datatypes=datatypes, ylab='fraction')
    fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
    plt.close(fig)
    # kh and kv
    kgroups = list(set(khv_rch_points.loc[khv_rch_points['dtype']=='k','Group']))
    for kidx  in ['kh', 'kv']:

        # ks layers in same plot (per group) (tricky)
        for g in kgroups:
            title = 'Log10 {} all layers {}'.format(kidx,g)
            print(title)
            sites = {}
            layers = []
            for l in range(smt.layers):
                sites['{} L{:02d}'.format(g,l)] = khv_rch_points.loc[khv_rch_points['Group']==g,'id']
                layers.append(l)
            datatypes = [kidx for e in sites.keys()]

            fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                                   layers=layers, supergroup=title,
                                                   sites=sites,
                                                   datatypes=datatypes, ylab='log({})'.format(kidx),
                                                   transformation=_log10)
            fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
            plt.close(fig)

        # ks groups in same plot (per layer)(tricky)
        for l in range(smt.layers):
            title = 'Log10 {} layer {:02d}'.format(kidx,l)
            print(title)
            sites = {}
            layers = []
            for g in kgroups:
                sites[g] = khv_rch_points.loc[khv_rch_points['Group']==g,'id']
                layers.append(l)
            datatypes = [kidx for e in sites.keys()]

            fig, ax = plot_supergroup_ppp_boxplots(filter_strs=filter_strs, param_nc=param_nc,
                                                   layers=layers, supergroup=title,
                                                   sites=sites,
                                                   datatypes=datatypes, ylab='log({})'.format(kidx),
                                                   transformation=_log10)
            fig.savefig(os.path.join(outdir, title.replace(' ', '_') + '.png'))
            plt.close(fig)





if __name__ == '__main__':
    plot_all_ppp_boxplots(r"T:\Temp\temp_gw_files\test_ppp_plots",'filter1')