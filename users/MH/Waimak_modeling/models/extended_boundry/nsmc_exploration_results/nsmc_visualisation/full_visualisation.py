# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/11/2017 6:30 PM
"""

from __future__ import division
from core import env
import os
import time
from observations.obs_2d_mean_sd import plot_all_2d_obs
from observations.obs_box_plots import plot_obs_all_all_boxplots
from observations.plot_ucn import plot_all_cons
from parameters.param_ppp_boxplots import plot_all_ppp_boxplots
from parameters.spatially_average_param import plt_all_spatial_param
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.end_member_mixing_analysis.visualise_filter import \
    plot_river_mixing
from traceback import format_exc


def full_vis(outdir, filters_strs):
    try:
        plot_all_2d_obs(os.path.join(outdir, 'obs_2d'), filters_strs)
    except:
        print(format_exc())

    try:
        plot_obs_all_all_boxplots(os.path.join(outdir, 'obs_1d'), filters_strs)
    except:
        print(format_exc())

    try:
        plot_all_cons(os.path.join(outdir, 'concentrations'), filters_strs)
    except:
        print(format_exc())

    try:
        plot_all_ppp_boxplots(os.path.join(outdir, 'param_1d'), filters_strs)
    except:
        print(format_exc())

    try:
        plt_all_spatial_param(os.path.join(outdir, 'param_2d'), filters_strs)
    except:
        print(format_exc())

    try:
        plot_river_mixing(os.path.join(outdir, 'emma_plots'), filters_strs)
    except:
        print(format_exc())


if __name__ == '__main__':
    filter_lists = [  # todo
        ['emma_converge', '~0emma_converge'],
        ['n_converge', '~0n_converge'],
        ['emma_ewf_wt', 'run_mt3d'],

        ['emma_chch_wt', 'run_mt3d'],
        ['emma_str_wt', 'run_mt3d'],
        ['emma_no_wt', 'run_mt3d'],
        ['emma_eq_wt', 'run_mt3d'],
        ['emma_ewf_wt', 'emma_chch_wt'],
        ['emma_no_wt', 'emma_eq_wt'],
        ['emma_chch_wt', 'emma_str_wt']

    ]

    base_dir = env.gw_met_data("mh_modeling/plots")
    t = time.time()
    for i, fstr in enumerate(filter_lists):
        print '{}: {} of {} {} min from start'.format(fstr, i, len(filter_lists), (t - time.time()) / 60)
        full_vis(os.path.join(base_dir, '_'.join(fstr).replace('~', 'not')), fstr)
