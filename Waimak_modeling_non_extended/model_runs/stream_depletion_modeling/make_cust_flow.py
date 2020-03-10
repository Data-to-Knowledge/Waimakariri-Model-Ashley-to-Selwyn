"""
Author: matth
Date Created: 9/06/2017 9:06 AM
"""

from __future__ import division
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
import os
from stream_depletion_model_setup import setup_and_run_stream_dep, setup_and_run_stream_dep_multip
import transient_inputs as tinputs
import numpy as np
from base_str_dep_runs import get_starting_heads_for_sd
from copy import deepcopy
import multiprocessing
import logging


def setup_runs():
    """
    function to return the model path for the fully naturalized fully transient run and ensure that the model has been
    run this run has a ss and then full year run from july to end of june  ## for SD 150
    :param recalc:
    :return: path to the model (without and extension)
    """
    #set up fully transient run
    base_path = "{}/add_flow_to_cust".format(base_mod_dir)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

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

    kwargs = {'name': None,
              'base_dir': base_path,
              'stress_vals': spv,
              'stress_to_month': stress_to_month,
              'wells_to_turn_on': {0: []},
              'ss': ss,
              'sy': sy,
              'silent': False,
              'solver': 'pcg',
              'model': 'mfnwt',  # solver will be changed internally to nwt
              'start_heads': start_hds,
              'flow_to_add_to_cust': 0}

    runs = []
    flows_to_add = [0.1,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,7,10]
    for flow in flows_to_add:
        temp = deepcopy(kwargs)
        temp['name'] = 'v_{}_added'.format(flow)
        temp['flow_to_add_to_cust'] = flow * 86400
        runs.append(temp)


    return runs

if __name__ == '__main__':
    runtype = 1
    if runtype == 1:
        multiprocessing.log_to_stderr(logging.DEBUG)
        runs = setup_runs()
        pool_size = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=pool_size)
        pool_outputs = pool.map(setup_and_run_stream_dep_multip, runs)
        pool.close()  # no more tasks
        pool.join()
    elif runtype == 2: # debug run
        runs = setup_runs()
        setup_and_run_stream_dep(**runs[2])
