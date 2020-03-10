"""
Author: matth
Date Created: 1/06/2017 11:16 AM
"""

from __future__ import division
import os
import numpy as np
from transient_run import setup_and_run_transient_multip
import transient_inputs as tinputs
import multiprocessing
import logging
from copy import deepcopy
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
import users.MH.Waimak_modeling.model_tools as mt
from plot_heads import plt_transient_heads
from glob import glob


def setup_all_runs():

    spv = {'nper': 13,
           'perlen': 30,
           'nstp': 2,
           'steady': [True, False, False, False, False, False, False, False, False, False, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: None,
                       1: 7,
                       2: 8,
                       3: 9,
                       4: 10,
                       5: 11,
                       6: 12,
                       7: 1,
                       8: 2,
                       9: 3,
                       10: 4,
                       11: 5,
                       12: 6}
    base_path = "{}/sd_comparative_runs".format(base_mod_dir)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    base_kwargs = {'name': None,
                   'base_dir': base_path,
                   'stress_vals': spv,
                   'stress_to_month': stress_to_month,
                   'wells_to_turn_off': {0: []},
                   'ss': ss,
                   'sy': sy,
                   'silent': True,
                   'solver': 'pcg',  # internally changed to nwt
                   'model': 'mfnwt',
                   'start_heads': None,
                   'model_version': 'm_strong_vert',
                   'well_data': None,
                   'rch_data': None}

    runs = []

    # fully naturalized run without races
    fn_no_race = deepcopy(base_kwargs)
    fn_no_race['name'] = 'fully_nat_without_races'
    input_wells = mt.get_original_well_data()
    input_wells = input_wells[input_wells['type'] == 'boundry_flux']
    fn_no_race['well_data'] = input_wells
    fn_no_race['rch_data'] = mt.get_nat_rch()
    fn_no_race['str_type'] = 'nat'
    runs.append(fn_no_race)

    # full noise abstraction (consented volume)
    full_noise = deepcopy(base_kwargs)
    full_noise['name'] = 'fully_consented_abstraction'
    full_noise['well_data'] = mt.get_model_well_full_consented()
    runs.append(full_noise)

    # current abstraction
    ca = deepcopy(base_kwargs)
    ca['name'] = 'current_abstraction'
    runs.append(ca)

    # full noise abstraction (consented volume) with true cmp
    full_noise = deepcopy(base_kwargs)
    full_noise['name'] = 'fully_consented_abstraction_true_cmp'
    full_noise['well_data'] = mt.get_model_well_full_consented()
    fn_no_race['rch_data'] = mt.get_true_cmp_rch()
    runs.append(full_noise)

    # current abstraction with true cmp
    ca = deepcopy(base_kwargs)
    ca['name'] = 'current_abstraction_true_cmp'
    fn_no_race['rch_data'] = mt.get_true_cmp_rch()
    runs.append(ca)

    # current without near field stream depleting wells
    # this is not run and may be depreciated
    #ca_no_near_field = deepcopy(base_kwargs)
    #ca_no_near_field['name'] = 'current_abstraction_without_near_field_sd_wells'
    #ca_no_near_field['wells_to_turn_off'][0] = []  #  need a list of near field wells
    #runs.append(ca_no_near_field)


    # GMP
    gmp = deepcopy(base_kwargs)
    gmp['name'] = 'gmp_simulation'
    gmp['rch_data'] = mt.get_gmp_rch()
    runs.append(gmp)

    # GMP + full noise abstraction
    gmp_full_ab = deepcopy(base_kwargs)
    gmp_full_ab['name'] = 'gmp_full_consented_abstraction'
    gmp_full_ab['rch_data'] = mt.get_gmp_rch()
    gmp_full_ab['well_data'] = mt.get_model_well_full_consented()
    runs.append(gmp_full_ab)

    return runs

def start_process():
    print('Starting', multiprocessing.current_process().name)

def ex_runs_via_mprocessing():
    multiprocessing.log_to_stderr(logging.DEBUG)
    runs = setup_all_runs()
    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    pool_outputs = pool.map(setup_and_run_transient_multip, runs)
    pool.close()  # no more tasks
    pool.join()
    with open(r"C:\Users\matth\Desktop\run_status.txt",'w') as f:
        wr = ['{}: {}\n'.format(e[0],e[1]) for e in pool_outputs]
        f.writelines(wr)

def plt_result_hds():

    out_plot_dir = "{}/sd_comparative_runs/plots".format(base_mod_dir)
    if not os.path.exists(out_plot_dir):
        os.makedirs(out_plot_dir)
    head_files = glob("{}/sd_comparative_runs/*/*.hds".format(base_mod_dir))
    stress_to_month = {0: None,
                       1: 7,
                       2: 8,
                       3: 9,
                       4: 10,
                       5: 11,
                       6: 12,
                       7: 1,
                       8: 2,
                       9: 3,
                       10: 4,
                       11: 5,
                       12: 6}
    for path in head_files:
        outdir = '{}/{}'.format(out_plot_dir,os.path.basename(path).split('.')[0])
        plt_transient_heads(path,outdir,[0,1,2,3],stress_to_month)




if __name__ == '__main__':
    ex_runs_via_mprocessing()
    plt_result_hds()