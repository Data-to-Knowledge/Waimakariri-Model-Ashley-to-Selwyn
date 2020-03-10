"""
Author: matth
Date Created: 30/05/2017 12:42 PM
"""

from __future__ import division
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, base_mod_dir2
import users.MH.Waimak_modeling.model_tools as mt
import os
import multiprocessing
import logging
from stream_depletion_model_setup import setup_and_run_stream_dep, setup_and_run_stream_dep_multip
from copy import deepcopy
import transient_inputs as tinputs
import numpy as np
import socket
from base_str_dep_runs import get_starting_heads_for_sd
import time
from future.builtins import input

def setup_runs(well_list,add_to_cust=False):
    start_hds = get_starting_heads_for_sd()
    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    spv = {'nper': 5,
           'perlen': 30,
           'nstp': 2,
           'steady': [False, False, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: 10,
                       1: 11,
                       2: 12,
                       3: 1,
                       4: 2}

    if add_to_cust:
        base_path = "{}/well_by_well_str_dep_sd150_1cu_added_to_cust".format(base_mod_dir2)
    else:
        base_path = "{}/well_by_well_str_dep_sd150".format(base_mod_dir)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    if os.path.exists('{}/error_log.txt'.format(base_path)):
        os.remove('{}/error_log.txt'.format(base_path))

    if add_to_cust:
        add_h20 = 86400
    else:
        add_h20 = 0
    base_kwargs = {'name': None,
                   'base_dir': base_path,
                   'stress_vals': spv,
                   'stress_to_month': stress_to_month,
                   'wells_to_turn_on': {0: []},
                   'ss': ss,
                   'sy': sy,
                   'silent': True,
                   'solver': 'pcg',
                   'start_heads': start_hds,
                   'model': 'mfnwt',
                   'flow_to_add_to_cust': add_h20}  # solver will be changed internally to nwt

    out_runs = []
    for well in well_list:
        temp_kwargs = deepcopy(base_kwargs)
        temp_kwargs['wells_to_turn_on'][0] = [well]
        temp_kwargs['name'] = 'turn_on_{}'.format(well.replace('/','_'))
        out_runs.append(temp_kwargs)

    return out_runs


def start_process():
    print('Starting', multiprocessing.current_process().name)


def well_by_well_depletion(well_list,add_to_cust=False):
    t = time.time()
    multiprocessing.log_to_stderr(logging.DEBUG)
    runs = setup_runs(well_list,add_to_cust)
    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    pool_outputs = pool.map(setup_and_run_stream_dep_multip, runs)
    pool.close()  # no more tasks
    pool.join()
    with open(r"{}/well_by_well run_status.txt".format(base_mod_dir), 'w') as f:
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in pool_outputs]
        f.writelines(wr)
    print('{} runs completed in {} minutes'.format(len(well_list),((time.time()-t)/60)))


if __name__ == '__main__':
    _comp = socket.gethostname()
    cont = input(
        'running well_by_well_str_dep:\n continue y/n\n')
    if cont != 'y':
        raise ValueError('script aborted so as not to overwrite')

    if _comp == 'HP1639':
        runs = setup_runs(mt.get_n_wai_well_list()[0:1],True)
        runs[0]['silent'] = False
        setup_and_run_stream_dep(**runs[0])
    elif _comp == 'DHI-Runs02':
        well_list = mt.get_n_wai_well_list()[0:601] #601 #this split assumes 16 cores on dhiruns 02 and 24 v cores on gwater02
        #well_by_well_depletion(well_list,False) already ran, don't need to re-run
        well_by_well_depletion(well_list, True)
    elif _comp == 'GWATER02':
        well_list = mt.get_n_wai_well_list()[601:]
        #well_by_well_depletion(well_list,False) already ran, don't need to re-run
        well_by_well_depletion(well_list, True)
    else:
        raise ValueError('unexpected computer'.format(_comp))
