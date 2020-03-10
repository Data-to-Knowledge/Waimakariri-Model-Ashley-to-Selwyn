# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/06/2018 9:02 AM
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy
import os


def plot_comparison_boxplot(data, gmp_data, outpath, scen_names, redlevels):
    """
    :param data: series with multiindex (scenario,stat)
    :param outpath:
    :return:
    """

    fig, ax = plt.subplots(figsize=(18.5, 9.5))
    positions = [1]
    fake_data = [[1, 2, 3]]
    keys = [None]
    labels = ['GMP']
    x = 2
    for red in redlevels:
        for scen in scen_names:
            positions.append(deepcopy(x))
            fake_data.append([1, 2, 3])
            keys.append('{}_red{}'.format(scen,red))
            labels.append(scen)
            x += 0.5
        x += 0.5
    ymax = []
    ymin = []
    box_plot = ax.boxplot(fake_data,
                          positions=positions,
                          labels=labels,
                          widths=0.25)
    for box_no, key in enumerate(keys):
        if box_no == 0:
            q1_start = gmp_data.loc['5%']
            q2_start = gmp_data.loc['25%']
            q3_start = gmp_data.loc['50%']
            q4_start = gmp_data.loc['75%']
            q4_end = gmp_data.loc['95%']
        else:
            q1_start = data.loc[key, '5%']
            q2_start = data.loc[key, '25%']
            q3_start = data.loc[key, '50%']
            q4_start = data.loc[key, '75%']
            q4_end = data.loc[key, '95%']

        ymax.append(q4_end * 1.1)
        ymin.append(q1_start * 0.9)

        # Lower cap
        box_plot['caps'][2 * box_no].set_ydata([q1_start, q1_start])
        # xdata is determined by the width of the box plot

        # Lower whiskers
        box_plot['whiskers'][2 * box_no].set_ydata([q1_start, q2_start])

        # Higher cap
        box_plot['caps'][2 * box_no + 1].set_ydata([q4_end, q4_end])

        # Higher whiskers
        box_plot['whiskers'][2 * box_no + 1].set_ydata([q4_start, q4_end])

        # Box
        box_plot['boxes'][box_no].set_ydata([q2_start,
                                             q2_start,
                                             q4_start,
                                             q4_start,
                                             q2_start])

        # Median
        box_plot['medians'][box_no].set_ydata([q3_start, q3_start])

    # The y axis is rescaled to fit the new box plot completely with 10%
    # of the maximum value at both ends
    ax.set_ylim([min(ymin), max(ymax)])
    ax.set_title(data.name)
    ax.set_ylabel('g/m3 N')
    an_y = .04
    ax.annotate('10% reduction',(.31,an_y), xycoords='figure fraction',weight = 'bold')
    ax.annotate('20% reduction',(.53,an_y), xycoords='figure fraction',weight = 'bold')
    ax.annotate('30% reduction',(.76,an_y), xycoords='figure fraction',weight = 'bold')
    fig.savefig(outpath)
    plt.close(fig)


if __name__ == '__main__':
    # do waimak and interzone plots.
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\reduction_plots\waimak_zone"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    alldata = pd.read_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
        index_col=0, header=[0, 1])
    allgmp = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
                         index_col=0, header=[0,1])['gmp']
    for site in list(alldata.index):
        print(site)
        plot_comparison_boxplot(data=alldata.loc[site],
                                gmp_data=allgmp.loc[site],
                                outpath=os.path.join(outdir,'{}.png'.format(site)),
                                scen_names=['dairy', 'kgha_15', 'kgha_5'],
                                redlevels=[10, 20, 30])

    print('done')

    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\reduction_plots\interzone"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    alldata = pd.read_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\interzone\all_n_interzone.csv",
        index_col=[0,1], header=[0, 1])
    allgmp = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\interzone\all_n_interzone.csv",
                         index_col=[0,1], header=[0,1])['gmp']
    for site in list(alldata.index):
        print(site)
        plot_comparison_boxplot(data=alldata.loc[site],
                                gmp_data=allgmp.loc[site],
                                outpath=os.path.join(outdir,'{}.png'.format('_'.join(site))),
                                scen_names=['dairy', 'kgha_15', 'kgha_5'],
                                redlevels=[10, 20, 30])

    print('done')
