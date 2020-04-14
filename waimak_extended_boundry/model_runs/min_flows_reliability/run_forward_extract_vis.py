# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/10/2017 4:29 PM
"""

from __future__ import division
import os
import time

from waimak_extended_boundry.model_run_tools.forward_quanity_support.forward_runs import setup_run_args, run_forward_runs
from waimak_extended_boundry.model_run_tools.forward_quanity_support.extract_data_for_forward_runs import gen_all_outdata_forward_runs, extract_and_save_all_cc_mult_missing_w
from visualise_data_from_fruns import plot_and_save_forward_vis

if __name__ == '__main__':
    # a convinence set to run everything for the forward runs
    # note that the well reliablity scripts where not completed at the time that this was run so are not included here
    # and must be run separately, if I need to re-run the scripts at a later date include this item

    # done on rdsprod03
    #### inputs to define for each run####
    safemode = True
    model_id = 'StrOpt'
    model_dir_path = r"D:\mh_model_runs\{}_forward_runs_2017_10_29".format(model_id)  # path on rdsprod03
    results_dir = r"D:\mh_model_runs\{}_forward_runs_2017_10_29_results".format(model_id)
    cc_to_waimak_only = True
    notes = """ 
    Naturalisation applied to full domain, CC component of LSR changes and pc5 applied to only waimakariri, 
    but ccmult and missing h20 is waimak only these senarios are without the n carpet drains unless 
    explecitly stated in the model name (_w_ncar)
    pumping changes only applied to Waimakariri
    run in {}
    """.format(model_dir_path)
    run_models = False # a quick way to only run part of this script if something breaks
    amalg_results = True # a quick way to only run part of this script if someting breaks

    #### run the models ####
    if run_models:
        if safemode:
            if os.path.exists(model_dir_path):
                cont = input(
                    'run all forward runs, this could overwrite item in :\n {} \n continue y/n\n'.format(
                        model_dir_path)).lower()
                if cont != 'y':
                    raise ValueError('script aborted so as not to potentially overwrite {}'.format(model_dir_path))

        runs = setup_run_args(model_id, model_dir_path, cc_to_waimak_only=cc_to_waimak_only)

        t = time.time()
        run_forward_runs(runs, model_dir_path, notes)
        print('{} runs in __ min'.format(len(runs)))
        print((time.time() - t) / 60)

    if amalg_results:
        gen_all_outdata_forward_runs(
            model_dir_path,
            results_dir,
            False)
        extract_and_save_all_cc_mult_missing_w(model_dir_path, os.path.join(results_dir, "ccmult_extract.csv"))
        print('done')

    plot_and_save_forward_vis(
        os.path.join(results_dir, "overview_plots"),
        os.path.join(results_dir, "{}_relative_data.csv".format(model_id)),
        os.path.join(results_dir, "{}_meta_data.csv".format(model_id))
    )
