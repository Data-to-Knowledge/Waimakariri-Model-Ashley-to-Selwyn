"""
Author: matth
Date Created: 7/09/2017 5:19 PM
"""

from __future__ import division
import os
import multiprocessing
import logging
from waimak_extended_boundry.model_run_tools.stream_depletion_support.stream_depletion_model_setup import setup_and_run_stream_dep_multip
from copy import copy
import time
from waimak_extended_boundry.model_run_tools.model_bc_data.starting_hds_ss_sy import get_ss_sy, get_starting_heads_sd30, get_sd_well_list
from waimak_extended_boundry import smt
import psutil
import datetime
from future.builtins import input
from waimak_extended_boundry.model_run_tools.stream_depletion_support.base_sd_runs import get_sd_spv


def setup_runs_sd30(model_id, well_list, base_path, ss, sy, start_heads):
    """
    sets up the model runs for stream depletion 30 day assessment
    :param model_id: the NSMC realisation to use
    :param well_list: the list of wells to assess stream depletion one model will be run for each well
    :param base_path: the path to put all of the folders containing each well model
    :param ss: specific storage for the model (either integer or k,i,j array)
    :param sy: specific yield for the model (either integer or k,i,j array)
    :param start_heads: the starting heads for the model (k,i,j array)
    :return:
    """

    spv = get_sd_spv('sd30')

    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    if os.path.exists('{}/error_log.txt'.format(base_path)):
        os.remove('{}/error_log.txt'.format(base_path))

    base_kwargs = {
        'model_id': model_id,
        'base_dir': base_path,
        'stress_vals': spv,
        'ss': ss,
        'sy': sy,
        'silent': True,
        'start_heads': start_heads,
        'sd_7_150': 'sd30'}

    out_runs = []
    for well in well_list:
        temp_kwargs = copy(base_kwargs)
        temp_kwargs['wells_to_turn_on'] = {0: [well]}
        temp_kwargs['name'] = 'turn_on_{}_sd30'.format(well.replace('/', '_'))
        out_runs.append(temp_kwargs)

    return out_runs


def start_process():
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def well_by_well_depletion_sd30(model_id, well_list, base_path, notes, ss_sy_version=1):
    """
    run the well by well depletion for the 30 day stream depletion
    :param model_id: the NSMC realisation to use
    :param well_list: the list of wells to assess stream depletion one model will be run for each well
    :param base_path: the path to put all of the folders containing each well model
    :return:
    """
    if os.path.exists(base_path):
        cont = input("the base path already exists: \n {}\n do you want to continue y/n\n".format(base_path))
        if cont.lower() != 'y':
            raise KeyboardInterrupt('run  stopped to prevent overwrite of {}'.format(base_path))

    t = time.time()
    multiprocessing.log_to_stderr(logging.DEBUG)
    ss, sy = get_ss_sy(ss_sy_version=ss_sy_version)
    start_heads = get_starting_heads_sd30(model_id)
    runs = setup_runs_sd30(model_id, well_list, base_path, ss, sy, start_heads)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map_async(setup_and_run_stream_dep_multip, runs)
    while not results.ready():
        print('{} runs left of {}'.format(results._number_left, len(runs)))
        time.sleep(60 * 5)  # sleep 5 min between printing
    pool_outputs = results.get()
    pool.close()  # no more tasks
    pool.join()
    now = datetime.datetime.now()
    with open("{}/forward_run_log/{}_SD30_run_status_{}_{:02d}_{:02d}_{:02d}_{:02d}.txt".format(smt.sdp, model_id,
                                                                                                now.year, now.month,
                                                                                                now.day, now.hour,
                                                                                                now.minute), 'w') as f:
        f.write(str(notes) + '\n')
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in pool_outputs]
        f.writelines(wr)
    print('{} runs completed in {} minutes'.format(len(well_list), ((time.time() - t) / 60)))


if __name__ == '__main__':
    notes = """ """
    model_id = 'opt'
    well_list = get_sd_well_list(model_id)
    base_path = r"C:\Users\MattH\Desktop\test_sd30"
    well_by_well_depletion_sd30(model_id, well_list, base_path, notes)
    print('done')
