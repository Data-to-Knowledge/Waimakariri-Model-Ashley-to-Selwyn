"""
Author: matth
Date Created: 29/10/2017 12:02 PM
"""

from __future__ import division
from core import env
import os
import flopy
from glob import glob
import multiprocessing
import logging
import psutil
import time
from users.MH.Waimak_modeling.supporting_data_path import sdp
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.convergance_check import modflow_converged
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools import zip_non_essential_files



def _runnwt_model_forward_runs(path):
    try:
        exe_name = "{}/models_exes/MODFLOW-NWT_1.1.2/MODFLOW-NWT_1.1.2/bin/MODFLOW-NWT_64.exe".format(sdp)
        success, buff = flopy.mbase.run_model(exe_name, os.path.basename(path), os.path.dirname(path))
        if success:
            con = modflow_converged(path.replace('.nam', '.list'))
            zip_non_essential_files(os.path.dirname(path), include_list=True)

    except Exception as val:
        print(path, val)

def start_process():
    """
    function to run at the start of each multiprocess sets the priority lower
    :return:
    """
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def rerun_modflow_mem_errors(base_dir):
    paths = glob(os.path.join(base_dir, '*/*.nam'))
    run_paths = []
    for path in paths:
        if os.path.exists(os.path.join(os.path.dirname(path),'non_essential_components.zip')):
            continue

        if not os.path.exists(path.replace('.nam','.hds')):
            continue

        if os.path.getsize(path.replace('.nam','.hds')) != 0:
            continue

        run_paths.append(path)

    multiprocessing.log_to_stderr(logging.DEBUG)
    pool_size = int(multiprocessing.cpu_count()/4*3)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map_async(_runnwt_model_forward_runs, run_paths)
    while not results.ready():
        print('{} runs left of {}'.format(results._number_left, len(run_paths)))
        time.sleep(60*5)  # sleep 5 min between printing
    pool_outputs = results.get()
    pool.close()  # no more tasks
    pool.join()
    print(pool_outputs)

if __name__ == '__main__':
    rerun_modflow_mem_errors(r"D:\mh_model_runs\StrOpt_forward_runs_2017_10_29")