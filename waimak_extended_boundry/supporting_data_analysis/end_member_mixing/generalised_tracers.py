# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 10/11/2017 8:55 AM
"""

from __future__ import division
import numpy as np
import pandas as pd
import itertools
import os
from scipy.stats import mode
import matplotlib.pyplot as plt
from future.builtins import input


def mc_calc_end_members(outdir, sites, means, sds, tracers, end_members, a_errors, n=10000, plot=True):
    if len(end_members) > 4:
        cont = input('this will probably cause a memory error continue y/n')
        if cont != 'y':
            raise ValueError('user stopped to prevent memory error')
    percentages = np.array(list(itertools.product(*[range(101) for e in end_members])))
    percentages = percentages[percentages.sum(axis=1) == 100] / 100
    cons = {}
    for trace in tracers:
        for i, end in enumerate(endmembers):
            if i == 0:
                temp = np.random.normal(means[trace][end], sds[trace][end], n)[np.newaxis, np.newaxis,
                       :] * percentages[:, i, np.newaxis]
            else:
                temp += np.random.normal(means[trace][end], sds[trace][end], n)[np.newaxis, np.newaxis,
                        :] * percentages[:, i, np.newaxis]
        cons[trace] = temp[0]

    # shape will be (percentages,component(coastal,inland,river), itterations)
    outdata_median = pd.DataFrame(index=sites.keys(), columns=end_members)
    outdata_5th = pd.DataFrame(index=sites.keys(), columns=end_members)
    outdata_95th = pd.DataFrame(index=sites.keys(), columns=end_members)
    outdata_mode = pd.DataFrame(index=sites.keys(), columns=end_members)
    outdata_mean = pd.DataFrame(index=sites.keys(), columns=end_members)
    plot_dir = os.path.join(outdir, 'plots')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    for site in sites.keys():
        idx = np.ones(cons.values()[0].shape).astype(bool)
        for trace in tracers:
            upper = sites[site]['{}_upper'.format(trace)] + a_errors[trace]
            lower = sites[site]['{}_lower'.format(trace)] - a_errors[trace]
            idx = idx & ((cons[trace] <= upper) & (cons[trace] >= lower))

        temp = []
        for i in range(idx.shape[1]):
            if not idx[:, i].any():
                continue
            temp.append(percentages[idx[:, i]])
        temp_out = np.concatenate(temp, axis=0)
        out_5th = np.percentile(temp_out, 5, axis=0)
        outdata_5th.loc[site] = out_5th
        out_50th = np.percentile(temp_out, 50, axis=0)
        outdata_median.loc[site] = out_50th
        out_95th = np.percentile(temp_out, 95, axis=0)
        outdata_95th.loc[site] = out_95th
        outdata_mode.loc[site] = mode(temp_out, axis=0).mode[0]
        outdata_mean.loc[site] = temp_out.mean(axis=0)
        if plot:
            fig, ax = plt.subplots(figsize=(18.5, 9.5))
            for i, end in enumerate(end_members):
                ax.hist(temp_out[:, i], bins=101, label=end, alpha=0.5)
            ax.set_title(site)
            ax.legend()
            fig.savefig(os.path.join(plot_dir, '{}.png'.format(site.replace('/', '_'))))
            plt.close(fig)

    outdata_5th.to_csv(os.path.join(outdir, '5th.csv'))
    outdata_median.to_csv(os.path.join(outdir, 'median.csv'))
    outdata_95th.to_csv(os.path.join(outdir, '95th.csv'))
    outdata_mode.to_csv(os.path.join(outdir, 'mode.csv'))
    outdata_mean.to_csv(os.path.join(outdir, 'mean.csv'))


if __name__ == '__main__':
    tracers = [
        'cl',
        'o18',
        'k'
    ]
    endmembers = [
        'inland',
        'coastal',
        'river',
        'eyre'
    ]

    means = {
        'cl': {
            'inland': 10.38,
            'coastal': 25.43,
            'river': 1.05,
            'eyre': 3.67
        },
        'o18': {
            'inland': -8.76,
            'coastal': -8.00,
            'river': -9.25,
            'eyre': -8.90
        },
        'k': {
            'inland': 1.23,
            'coastal': 1.67,
            'river': 0.57,
            'eyre': 1.17
        }
    }

    sds = {
        'cl': {
            'inland': 4.21,
            'coastal': 9.49,
            'river': 0.15,
            'eyre': 0.29
        },
        'o18': {
            'inland': 0.21,
            'coastal': 0.64,
            'river': 0.22,
            'eyre': 0.54
        },
        'k': {
            'inland': 0.31,
            'coastal': 0.30,
            'river': 0.15,
            'eyre': 0.12
        }
    }

    a_errors = {
        'cl': 0.5,
        'o18': 0.2,
        'k': 0.15
    }
    # these are for most sites
    #targets = pd.read_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells\AdditionalTargets_5Tracers_Python.csv", index_col=0)
    # these are for a few more streams
    targets = pd.read_excel(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\OkokaKaiapoi_3Tracers_EMMAInputs.xlsx", index_col=0)
    sites = {}
    for site in targets.index:
        target = {}
        for trace in tracers:
            target['{}_lower'.format(trace)] = targets.loc[site, '{}_mean'.format(trace)] - targets.loc[site, '{}_stdev'.format(trace)]
            target['{}_upper'.format(trace)] = targets.loc[site, '{}_mean'.format(trace)] + targets.loc[site, '{}_stdev'.format(trace)]
        sites[site] = target

    mc_calc_end_members(
        outdir=r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells\4_members_tracers_k_18o_cl_new_streams",
        sites=sites,
        means=means,
        sds=sds,
        tracers=tracers,
        end_members=endmembers,
        a_errors=a_errors,
        plot=True,
        n=5000
    )
    print'done'
