"""
Author: matth
Date Created: 23/03/2018 4:16 PM
"""
from __future__ import division
import os
import time
from future.builtins import input
from forward_runs import setup_run_args, run_forward_runs
from waimak_extended_boundry.model_run_tools.forward_quanity_support.extract_data_for_forward_runs import gen_all_outdata_forward_runs, extract_and_save_all_cc_mult_missing_w
from visualise_data_from_fruns import plot_and_save_forward_vis
from waimak_extended_boundry.model_run_tools import trans_readme, get_stocastic_set
from waimak_extended_boundry.model_run_tools.model_setup.new_model_run import minium_new_model_run
from glob import glob
from traceback import format_exc
import socket

if __name__ == '__main__':

    # a convinence set to run everything for the forward runs
    # note that the well reliablity scripts where not completed at the time that this was run so are not included here
    # and must be run separately, if I need to re-run the scripts at a later date include this item
    # this is meant to easily run the forward runs for the non cc senarios

    # done on rdsprod03
    #### inputs to define for each run####
    safemode = False
    model_ids = get_stocastic_set()
    if socket.gethostname() == 'GWATER02':
        model_ids = get_stocastic_set()
        model_ids.reverse()
    run_modelses = [True for e in model_ids]
    date = '2018-03-24'
    base_dir = r"D:\mh_waimak_models\stocastic_forward"
    remove_carpet = True

    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    with open(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\forward_run_log\stocastic_forward_run_status_{}.txt".format(
                date), 'w') as f:
        pass

    #### run the models ####
    for i, (model_id, run_models) in enumerate(zip(model_ids, run_modelses)):
        print ("""
        ##########################################################
        starting forward runs for {}, the models will be re-run: {} model {} of {}
        ##########################################################
        """.format(model_id, run_models, i + 1, len(model_ids)))
        try:
            minium_new_model_run(model_id)
            model_dir_path = r"{}\{}_non_cc_forward_runs_{}".format(base_dir, model_id, date)
            results_dir = r"{}\{}_non_cc_forward_runs_{}_results".format(base_dir, model_id, date)
            cc_to_waimak_only = False
            if remove_carpet:
                add = 'carpet drains removed if not explicit in the model name (_w_ncar)'
            else:
                add = 'carpet drains not removed'
            notes = """ 
            Naturalisation changes applied to full domain, No climate change senarios run {}
            pumping changes only applied to Waimakariri, run in {}
            """.format(add, model_dir_path)
            if run_models:
                if safemode:
                    if os.path.exists(model_dir_path):
                        cont = input(
                            'run all forward runs, this could overwrite item in :\n {} \n continue y/n\n'.format(
                                model_dir_path)).lower()
                        if cont != 'y':
                            raise ValueError(
                                'script aborted so as not to potentially overwrite {}'.format(model_dir_path))

                if not os.path.exists(model_dir_path):
                    os.makedirs(model_dir_path)
                if not os.path.exists(results_dir):
                    os.makedirs(results_dir)
                runs = setup_run_args(model_id, model_dir_path, cc_to_waimak_only=cc_to_waimak_only, cc_runs=False)
                if not remove_carpet:
                    [e.update({'rm_ncarpet': False}) for e in runs]
                t = time.time()
                run_forward_runs(runs, model_dir_path, notes)
                print('{} runs in __ min'.format(len(runs)))
                print((time.time() - t) / 60)

            gen_all_outdata_forward_runs(
                model_dir_path,
                results_dir,
                False)
            extract_and_save_all_cc_mult_missing_w(model_dir_path, os.path.join(results_dir, "ccmult_extract.csv"))
            trans_readme(model_dir_path, results_dir)
            plot_and_save_forward_vis(
                os.path.join(results_dir, "overview_plots"),
                os.path.join(results_dir, "{}_relative_data.csv".format(model_id)),
                os.path.join(results_dir, "{}_meta_data.csv".format(model_id)),
                cc_runs=False, pc5_comp=True
            )

            # remove pickles from p drive
            delete_paths = glob(
                r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\pickled_files\temp_pickle_dir\model_{}_*.p".format(
                    model_id))
            for path in delete_paths:
                os.remove(path)
        except Exception as val:
            err = ('{}: {}\n'.format(model_id, format_exc().replace('\n', '')))
            with open(
                    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\forward_run_log\stocastic_forward_run_status_{}.txt".format(
                            date), 'a') as f:
                f.write(err)
