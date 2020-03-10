"""
Author: matth
Date Created: 22/05/2017 10:59 AM
"""

from __future__ import division
import numpy as np
from transient_run import setup_and_run_transient_multip
import transient_inputs as tiputs
import multiprocessing
import logging
from copy import copy
import os
import time


def setup_runs():
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
                       12: 6

                       }
    base_path = "D:/mattH/python_wm_runs/s_explore_transient"
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    base_kwargs = {'name': None,
                   'base_dir': base_path,
                   'stress_vals': spv,
                   'stress_to_month': stress_to_month,
                   'wells_to_turn_off': {0: []},
                   'ss': None,
                   'sy': None,
                   'silent': True,
                   'solver': 'pcg',
                   'model':'mfnwt'} # solver will be changed internally to nwt

    ss_vals = tiputs.get_ss_dist(100)
    sy_vals = tiputs.get_sy_dist(100)

    runs = []
    for i in range(100):
        temp = copy(base_kwargs)
        temp['ss'] = ss_vals[i]
        temp['sy'] = sy_vals[i]
        temp['name'] = 'run_{}_ss_{}_sy_{}'.format(i,ss_vals[i],sy_vals[i])
        runs.append(temp)
    # add investigation runs:
    ss_vals_2 = tiputs.get_ss_dist()
    sy_vals_2 = tiputs.get_sy_dist()
    for i in [1,5,10,25,50,75,90,95,99]:
        temp = copy(base_kwargs)
        temp['ss'] = np.percentile(ss_vals_2, i)
        temp['sy'] = np.percentile(sy_vals_2, i)
        temp['name'] = 'perc_{}_ss_{}_sy_{}'.format(i,ss_vals[i],sy_vals[i])
        runs.append(temp)

    return runs
def start_process():
    print('Starting', multiprocessing.current_process().name)

def ex_runs_via_mprocessing():
    multiprocessing.log_to_stderr(logging.DEBUG)
    runs = setup_runs()
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

if __name__ == '__main__':
    t = time.time()
    ex_runs_via_mprocessing()
    with open(r"C:\Users\matth\Desktop\runtime.txt",'w') as f:
        f.write('109 simulations took {} hours to run'.format((time.time() - t)/60/60))