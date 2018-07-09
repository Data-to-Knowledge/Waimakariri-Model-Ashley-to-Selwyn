# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/07/2018 11:53 AM
"""

from __future__ import division
from core import env
import numpy as np
from copy import deepcopy
import pandas as pd
import itertools
from percentage_reduction_maps import gen_stream_targets, gen_well_targets, get_current_pathway_n, gen_interzone_targets
from number_of_steps import lclasses, get_average_n, get_fractions, get_step_targets
import matplotlib.pyplot as plt
import os
import operator


def con_at_steps(step_reductions, pa_00, num_steps=10):
    """
    calculate the number of steps to reach a selected target scheme
    :param step_reductions: a dictionary of the percentage reduction to apply at each step e.g.
                            {'dairy_lowdrain' :20,
                            'sbda_lowdrain':10,
                            'dairy_highdrain':20,
                            'sbda_highdrain':5,
                            'lifestyle':0}
    :param pa_00: boolean if True, set teh pa rules to 0,0
    :param num_steps: number of steps to project
    :return:
    """
    assert isinstance(step_reductions, dict)
    assert set(step_reductions.keys()) == lclasses
    use_step_reductions = {}

    for k, v in step_reductions.items():
        use_step_reductions[k] = 1 - v / 100

    average_n = get_average_n()
    mode = {k: '50%' for k in average_n.keys()}
    current_paths = get_current_pathway_n(mode=mode, conservative_zones='use_mix', from_mt3d_runs=True,
                                          mt3d_add_pa=not pa_00, inc_interzone=True)

    fractions = get_fractions()
    outdata = pd.DataFrame(index=pd.MultiIndex.from_product((current_paths.keys(), range(0, num_steps + 1)),
                                                            names=['site', 'step']),
                           columns=[e + '_reduction' for e in lclasses] + ['n_con']
                           )
    for site in current_paths.keys():
        # calculate the average n concentration for the missing components (check this)
        n_proj = deepcopy(current_paths[site])
        n_proj_temp = deepcopy(current_paths[site])
        for lclass in lclasses:
            n_proj_temp += - average_n[site][lclass] * fractions[site][lclass]

        other_n_con = n_proj_temp / fractions[site]['other']

        # apply reduction, calculate new concentration
        for step in range(0, num_steps + 1):
            temp_use_reductions = {}
            for lclass in lclasses:  # assign the reduction to use and save it
                temp = step_reductions[lclass] * step / 100
                if temp >= 1:
                    temp = 1
                temp_use_reductions[lclass] = temp
                outdata.loc[(site, step), '{}_reduction'.format(lclass)] = temp
            new_n = other_n_con * fractions[site]['other']
            for lclass in lclasses:
                new_n += average_n[site][lclass] * (1 - temp_use_reductions[lclass]) * fractions[site][lclass]
            outdata.loc[(site, step), 'n_con'] = new_n

    return outdata


def run_scenarios(num_steps=10):
    pas = [True, False]

    # make reduction options
    dairy_options = [10, 25]
    sheep_options = [0, 5]
    lifestyle_options = [0]
    include_ld = [True, False]

    reduction_options = {}
    for d, s, ly, ld in itertools.product(dairy_options, sheep_options, lifestyle_options, include_ld):
        ld_name = 'exc'
        if ld:
            ld_name = 'inc'

        name = 'd{}_s{}_L{}_{}ld'.format(d, s, ly, ld_name)
        temp = {}
        temp['dairy_highdrain'] = d
        temp['sbda_highdrain'] = s
        temp['lifestyle'] = ly

        if ld:
            temp['dairy_lowdrain'] = d
            temp['sbda_lowdrain'] = s
        else:
            temp['dairy_lowdrain'] = 0
            temp['sbda_lowdrain'] = 0

        reduction_options[name] = temp
    scens = {}
    for red, pc5pa00 in itertools.product(reduction_options.keys(), pas):
        pa_name = 'without'
        if pc5pa00:
            pa_name = 'with'
        outname = '{}_{}_pc5pa00'.format(red, pa_name)
        scens[outname] = {'step_reductions': reduction_options[red],
                          'pa_00': pc5pa00}

    average_n = get_average_n()
    outdata = pd.DataFrame(index=pd.MultiIndex.from_product((average_n.keys(), range(0, num_steps + 1)),
                                                            names=['site', 'step']),
                           columns=pd.MultiIndex.from_product((scens.keys(),
                                                               [e + '_reduction' for e in lclasses] + ['n_con']),
                                                              names=['scenario', 'data'])
                           )
    for name, scen in scens.items():
        temp = con_at_steps(num_steps=num_steps, **scen)
        outdata[name] = temp
    return outdata


def plot_scen(data, site, keys_to_include, title=None, include_reduction_plots=True):
    """

    :param data: the data produced by run_scenarios
    :param site: which site to plot
    :param keys_to_include: keys: data key with PA removed e.g. d25_s0_L0_excld for d25_s0_L0_excld_with_pc5pa00
                            items: {'color':'b',
                                    'linestyle': '--',
                                    'label':'a label',
                                    kwargs for plot}
    :param title: the title of the plot
    :return: fig, axes
    """
    assert isinstance(keys_to_include, dict), 'keys_to_include must be dict... see docstring'
    if include_reduction_plots:
        fig, axes = plt.subplots(4, figsize=(8.27, 11.69), sharex=True)
    else:
        fig, axes = plt.subplots(2, figsize=(8.27, 5), sharex=True)
    pref_tar = get_step_targets('preferred')[site]
    alt_tar = get_step_targets('alternate')[site]
    maxes = [pref_tar, alt_tar]

    for key, kwargs in keys_to_include.items():
        assert isinstance(kwargs, dict), 'values of keys_to_include msut be dict, see doc'
        nopa_data = data['{}_with_pc5pa00'.format(key)].loc[site]
        pa_data = data['{}_without_pc5pa00'.format(key)].loc[site]
        maxes.append(pa_data.n_con.max())
        steps = pa_data.index.values

        # ax1 PC5 PA
        axes[0].plot(steps, pa_data.n_con, marker='.', **kwargs)

        # ax2 No PA
        axes[1].plot(steps, nopa_data.n_con, marker='.', **kwargs)

        if include_reduction_plots:
            # ax3 dairy reduction
            axes[2].plot(steps, nopa_data.dairy_highdrain_reduction * 100, marker='.', **kwargs)

            # ax4 spda reduction
            axes[3].plot(steps, nopa_data.sbda_highdrain_reduction * 100, marker='.', **kwargs)

    # make plot pretty
    # xaxis
    axes[-1].set_xlabel('reduction steps')
    axes[-1].set_xticks(steps)
    if max(maxes) > 5:
        spacer = 1
    elif max(maxes) > 2:
        spacer = 0.5
    else:
        spacer = 0.2
    # yaxis
    axes[0].set_ylabel('[N] with PA rules')
    axes[0].set_ylim([0, max(maxes) * 1.2])
    axes[0].set_yticks(np.arange(0, max(maxes)*1.2, spacer))

    axes[1].set_ylabel('[N] without PA rules')
    axes[1].set_ylim([0, max(maxes) * 1.2])
    axes[1].set_yticks(np.arange(0, max(maxes)*1.2, spacer))

    if include_reduction_plots:
        axes[2].set_ylabel('Dairy & dairy support reduction')
        axes[2].set_ylim([-5, 110])

        axes[3].set_ylabel('Other agricultural reduction')
        axes[3].set_ylim([-5, 110])

    handles, labels = axes[-1].get_legend_handles_labels()

    # sort legend them by labels

    hl = sorted(zip(handles, labels),
                key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)

    axes[-1].legend(handles2, labels2)

    # add limit lines (prefered and other)

    axes[0].axhline(pref_tar, c='r', ls=':', alpha=0.5)
    axes[0].text(.75, pref_tar, 'preferred target')
    axes[1].axhline(pref_tar, c='r', ls=':', alpha=0.5)
    axes[1].text(.75, pref_tar, 'preferred target')

    if pref_tar != alt_tar:
        axes[0].axhline(alt_tar, c='r', ls=':', alpha=0.5)
        axes[0].text(.75, alt_tar, 'COMAR target')
        axes[1].axhline(alt_tar, c='r', ls=':', alpha=0.5)
        axes[1].text(.75, alt_tar, 'COMAR target')


    if title is not None:
        axes[0].set_title(title)
    fig.tight_layout()

    return fig, axes


def plot_sites(data, outdir, sites=None, keys_to_include=None, include_reduciotn_plots=True):
    """
    plot the data for a select number of sites
    :param data: data from con_at_steps
    :param outdir: directory to save the data, created if it does not exist
    :param sites: sites to plot, if None then plot all sites
    :param keys_to_include: scenarios to include (see plot_scen for example), a default set is created when None is passed
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if sites is None:
        sites = get_average_n().keys()

    if keys_to_include is None:
        keys_to_include = {
            'd25_s0_L0_excld': {'color': 'k',
                                'linestyle': '--',
                                'label': 'D&DS: 25, SBDA: 0, ex. LD'},

            'd25_s0_L0_incld': {'color': 'k',
                                'linestyle': '-',
                                'label': 'D&DS: 25, SBDA: 0, inc. LD'},

            'd10_s5_L0_excld': {'color': 'b',
                                'linestyle': '--',
                                'label': 'D&DS: 10, SBDA: 5, ex. LD'},

            'd10_s5_L0_incld': {'color': 'b',
                                'linestyle': '-',
                                'label': 'D&DS: 10, SBDA: 5, inc. LD'},

        }

    for site in sites:
        if 'wdc' in site:
            title = 'WDC supply well at {}'.format(site.replace('wdc_','')).title()
        elif '_s' in site and '_shal' not in site:
            if 'kai' in site and 'court' not in site:
                title = '{} stream at {} rd.'.format(site.split('_')[0],site.split('_')[1]).title()
            else:
                title = '{} stream'.format(site.split('_')[0]).title()
        elif 'interzone' in site.lower():
            title = site.replace('_', ' ').title()
        else:
            title = 'Private wells near {}'.format(site.replace('_', ', ')).title()

        fig, axes = plot_scen(data, site, keys_to_include, title=title,
                              include_reduction_plots=include_reduciotn_plots)
        fig.savefig(os.path.join(outdir, '{}.png'.format(site)))
        plt.close(fig)


if __name__ == '__main__':
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\step_projections"
    test = run_scenarios()
    plot_sites(test, outdir=outdir)
    plot_sites(test, outdir=os.path.join(outdir,'slimmed_plots'),include_reduciotn_plots=False)
    test.to_csv(os.path.join(outdir,'all_step_data.csv'))
    print('done')
