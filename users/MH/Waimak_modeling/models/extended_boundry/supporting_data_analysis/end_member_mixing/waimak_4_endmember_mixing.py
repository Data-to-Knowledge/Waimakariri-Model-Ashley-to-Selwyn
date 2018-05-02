# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/11/2017 4:42 PM
"""
# -*- coding: utf-8 -*-

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import itertools
import os
from scipy.stats import mode
import matplotlib.pyplot as plt


def mc_calc_end_members(outdir, sites, o18_u, o18_s, cl_u, cl_s, n=10000, export_raw=False):
    percentages = np.array(list(itertools.product(range(101), range(101), range(101), range(101))))
    percentages = percentages[percentages.sum(axis=1) == 100] / 100

    o18 = np.random.normal(o18_u['inland'], o18_s['inland'], n)[np.newaxis, np.newaxis, :]* percentages[:,0, np.newaxis]
    o18 += np.random.normal(o18_u['coastal'], o18_s['coastal'], n)[np.newaxis, np.newaxis, :] * percentages[:,1, np.newaxis]
    o18 += np.random.normal(o18_u['river'], o18_s['river'], n)[np.newaxis, np.newaxis, :] * percentages[:,2, np.newaxis]
    o18 += np.random.normal(o18_u['eyre'], o18_s['eyre'], n)[np.newaxis, np.newaxis, :] * percentages[:,3, np.newaxis]
    o18 = o18[0]

    cl = np.random.normal(cl_u['inland'], cl_s['inland'], n)[np.newaxis, np.newaxis, :] * percentages[:, 0, np.newaxis]
    cl += np.random.normal(cl_u['coastal'], cl_s['coastal'], n)[np.newaxis, np.newaxis, :] * percentages[:, 1, np.newaxis]
    cl += np.random.normal(cl_u['river'], cl_s['river'], n)[np.newaxis, np.newaxis, :] * percentages[:, 2, np.newaxis]
    cl += np.random.normal(cl_u['eyre'], cl_s['eyre'], n)[np.newaxis, np.newaxis, :] * percentages[:, 3, np.newaxis]
    cl = cl[0]

    # shape will be (percentages,component(coastal,inland,river), itterations)
    outdata_median = pd.DataFrame(index=sites.keys(), columns=['coastal', 'inland', 'river', 'eyre'])
    outdata_5th = pd.DataFrame(index=sites.keys(), columns=['coastal', 'inland', 'river', 'eyre'])
    outdata_95th = pd.DataFrame(index=sites.keys(), columns=['coastal', 'inland', 'river', 'eyre'])
    outdata_mode = pd.DataFrame(index=sites.keys(), columns=['coastal', 'inland', 'river', 'eyre'])
    outdata_mean = pd.DataFrame(index=sites.keys(), columns=['coastal', 'inland', 'river', 'eyre'])
    plot_dir = os.path.join(outdir, 'plots')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    for site in sites.keys():
        cl_lower = sites[site]['cl_lower'] - 0.5  # these values are for the analytical error
        cl_upper = sites[site]['cl_upper'] + 0.5
        o18_lower = sites[site]['o18_lower'] - 0.2
        o18_upper = sites[site]['o18_upper'] + 0.2  # these values are for the analytical error
        idx = ((o18 <= o18_upper) &
               (o18 >= o18_lower) &
               (cl <= cl_upper) &
               (cl >= cl_lower))
        temp = []
        for i in range(idx.shape[1]):
            if not idx[:, i].any():
                continue
            temp.append(percentages[idx[:, i]])
        temp_out = np.concatenate(temp, axis=0)
        if export_raw: # save the original data
            export_dir = os.path.join(outdir,'raw_data')
            if not os.path.exists(export_dir):
                os.makedirs(export_dir)
            saver = pd.DataFrame(temp_out, columns=['coastal', 'inland', 'river', 'eyre'])
            for key in ['coastal', 'inland', 'river', 'eyre']:
                saver.loc[key].to_csv(os.path.join(export_dir,'{}_{}.txt'.format(site,key)),
                                        index=False, header=False, sep=' ')

        out_5th = np.percentile(temp_out, 5, axis=0)
        outdata_5th.loc[site] = out_5th
        out_50th = np.percentile(temp_out, 50, axis=0)
        outdata_median.loc[site] = out_50th
        out_95th = np.percentile(temp_out, 95, axis=0)
        outdata_95th.loc[site] = out_95th
        outdata_mode.loc[site] = mode(temp_out, axis=0).mode[0]
        outdata_mean.loc[site] = temp_out.mean(axis=0)

        fig, ax = plt.subplots(figsize=(18.5, 9.5))
        ax.hist(temp_out[:, 0], bins=101, color='orange', label='coastal', alpha=0.5)
        ax.hist(temp_out[:, 1], bins=101, color='r', label='inland', alpha=0.5)
        ax.hist(temp_out[:, 2], bins=101, color='b', label='river', alpha=0.5)
        ax.hist(temp_out[:, 3], bins=101, color='g', label='eyre', alpha=0.5)
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
    end_mean_o18 = {'inland': -8.76, 'coastal': -8.00, 'river': -9.25, 'eyre': -8.90}
    end_sd_o18 = {'inland': 0.21, 'coastal': 0.64, 'river': 0.22, 'eyre': 0.54}

    end_mean_cl = {'inland': 10.38, 'coastal': 25.43, 'river': 1.05, 'eyre': 3.67}
    end_sd_cl = {'inland': 4.21, 'coastal': 9.49, 'river': 0.15, 'eyre': 0.29}

    targets = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells\AdditionalTargets.csv",
        index_col=0)
    sites = {}
    for site in targets.index:
        sites[site] = {'o18_lower': targets.loc[site, 'o18_mean'] - targets.loc[site, 'o18_stdev'],
                       'o18_upper': targets.loc[site, 'o18_mean'] + targets.loc[site, 'o18_stdev'],
                       'cl_lower': targets.loc[site, 'cl_mean'] - targets.loc[site, 'cl_stdev'],
                       'cl_upper': targets.loc[site, 'cl_mean'] + targets.loc[site, 'cl_stdev']}
    mc_calc_end_members(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells\4_members",
        sites, end_mean_o18, end_sd_o18, end_mean_cl, end_sd_cl,n=5000)
    print'done'
