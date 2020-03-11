# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/11/2017 6:30 PM
"""

from __future__ import division
import os
import time
from observations.obs_2d_mean_sd import plot_all_2d_obs
from observations.obs_box_plots import plot_obs_all_all_boxplots
from observations.plot_ucn import plot_all_cons
from parameters.param_ppp_boxplots import plot_all_ppp_boxplots
from parameters.spatially_average_param import plt_all_spatial_param
from waimak_extended_boundry.nsmc_exploration_results.end_member_mixing_analysis.visualise_filter import \
    plot_river_mixing
from traceback import format_exc
import gc

def full_vis(outdir, filters_strs):
    try:
        plot_all_2d_obs(os.path.join(outdir, 'obs_2d'), filters_strs)
    except:
        print(format_exc())
    gc.collect()
    try:
        plot_obs_all_all_boxplots(os.path.join(outdir, 'obs_1d'), filters_strs)
    except:
        print(format_exc())
    gc.collect()
    try:
        plot_all_cons(os.path.join(outdir, 'concentrations'), filters_strs)
    except:
        print(format_exc())
    gc.collect()
    try:
        plot_all_ppp_boxplots(os.path.join(outdir, 'param_1d'), filters_strs)
    except:
        print(format_exc())
    gc.collect()
    try:
        plt_all_spatial_param(os.path.join(outdir, 'param_2d'), filters_strs)
    except:
        print(format_exc())
    gc.collect()
    try:
        plot_river_mixing(os.path.join(outdir, 'emma_plots'), filters_strs)
    except:
        print(format_exc())
    gc.collect()

if __name__ == '__main__':
    filter_lists = [
        ['emma_ewf_wt', 'run_mt3d'],

        ['emma_chch_wt', 'run_mt3d'],
        ['emma_str_wt', 'run_mt3d'],
        ['emma_no_wt', 'run_mt3d'],
        ['emma_eq_wt', 'run_mt3d'],
        ['emma_ewf_wt', 'emma_chch_wt'],
        ['emma_no_wt', 'emma_eq_wt'],
        ['emma_chch_wt', 'emma_str_wt'],
        ['emma_converge', '~0_emma_converge'],
        ['n_converge', '~0_n_converge']

    ]

    base_dir = "D:/mh_waimak_models/nsmc_plots"
    t = time.time()
    for i, fstr in enumerate(filter_lists):
        print '{}: {} of {} {} min from start'.format(fstr, i, len(filter_lists), (t - time.time()) / 60)
        full_vis(os.path.join(base_dir, '_'.join(fstr).replace('~', 'not')), fstr)
    print('finished {} filter sets after {} min'.format(len(filter_lists), (time.time() - t) / 60))
