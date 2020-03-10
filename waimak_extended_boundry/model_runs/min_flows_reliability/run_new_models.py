"""
Author: matth
Date Created: 24/03/2018 7:58 AM
"""

from __future__ import division
import multiprocessing
from waimak_extended_boundry import \
    get_model_name_path, get_stocastic_set
import os
import psutil
import logging


def start_process():
    """
    function to run at the start of each multiprocess sets the priority lower
    :return:
    """
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)

def _quick_function(model_id):
    print(model_id)
    try:
        get_model_name_path(model_id)
    except:
        pass

if __name__ == '__main__':
    model_ids = get_stocastic_set()[15:]

    multiprocessing.log_to_stderr(logging.DEBUG)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map(_quick_function, model_ids)
