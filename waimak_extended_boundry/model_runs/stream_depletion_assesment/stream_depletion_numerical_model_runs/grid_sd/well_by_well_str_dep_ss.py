# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 16/10/2017 1:45 PM
"""

from __future__ import division
from waimak_extended_boundry.model_run_tools.stream_depletion_support.grid_sd_support.ss_grid_sd_setup import grid_wells, setup_and_run_ss_grid_stream_dep_multip
import pandas as pd
import os
from copy import copy
import multiprocessing
import psutil
import time
import logging
import datetime
from waimak_extended_boundry.model_run_tools.model_bc_data.starting_hds_ss_sy import \
    get_sd_starting_hds
from waimak_extended_boundry import smt
from warnings import warn
from future.builtins import input

def setup_runs_grid(model_id, flux, base_path, start_heads, positive_flux_behavior='raise'):
    """
    sets up the model runs for stream depletion 7 day assessment
    :param model_id: the NSMC realisation to use
    :param flux: the flux to attribute to each well
    :param base_path: the path to put all of the folders containing each well model
    :param start_heads: the starting heads for the model (k,i,j array)
    :param positive_flux_behavior: 'raise' or 'warn'.  What to do on positive flux
    :return:
    """
    if positive_flux_behavior not in ['warn', 'raise']:
        raise ValueError('unexpected values for positive_flux_behviour: {}'.format(positive_flux_behavior))
    if flux > 0:
        if positive_flux_behavior == 'warn':
            warn('using a positive flux this means injection wells')
        elif positive_flux_behavior == 'raise':
            raise ValueError('flux is >= 0, i.e. injection, not permitted under current positive flux behavior')
    elif flux == 0:
        raise ValueError('why are you setting flux to zero!')

    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    base_kwargs = {
        'model_id': model_id,
        'base_dir': base_path,
        'silent': True,
        'start_heads': start_heads}

    out_runs = []
    well_list = grid_wells(flux)
    for well in well_list.index:
        temp_kwargs = copy(base_kwargs)
        temp_kwargs['wells_to_turn_on'] = pd.DataFrame(well_list.loc[well]).transpose()
        temp_kwargs['name'] = well
        out_runs.append(temp_kwargs)

    return out_runs


def start_process():
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def well_by_well_depletion_grid(model_id, flux, base_path, notes):
    """
    run the well by well depletion for the 7 day stream depletion
    :param model_id: the NSMC realisation to use
    :param flux: the flux to use for each well
    :param base_path: the path to put all of the folders containing each well model
    :param notes: a string to write to a text file
    :return:
    """
    if os.path.exists(base_path):
        cont = input("the base path already exists: \n {}\n do you want to continue y/n\n".format(base_path))
        if cont.lower() != 'y':
            raise KeyboardInterrupt('run  stopped to prevent overwrite of {}'.format(base_path))

    start_heads = get_sd_starting_hds(model_id, 'grid')
    t = time.time()
    multiprocessing.log_to_stderr(logging.DEBUG)
    runs = setup_runs_grid(model_id, flux, base_path, start_heads)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map_async(setup_and_run_ss_grid_stream_dep_multip, runs)
    while not results.ready():
        print('{} runs left of {}'.format(results._number_left, len(runs)))
        time.sleep(60 * 5)  # sleep 5 min between printing
    pool_outputs = results.get()
    pool.close()  # no more tasks
    pool.join()
    now = datetime.datetime.now()
    with open("{}/forward_run_log/{}_SDgrid_run_status_{}_{:02d}_{:02d}_{:02d}_{:02d}.txt".format(smt.sdp, model_id,
                                                                                                  now.year,
                                                                                               now.month, now.day,
                                                                                               now.hour, now.minute),
              'w') as f:
        f.write(str(notes) + '\n')
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in pool_outputs]
        f.writelines(wr)
    print('{} runs completed in {} minutes'.format(len(runs), ((time.time() - t) / 60)))


if __name__ == '__main__':
    setup_runs_grid('StrOpt',-10,'test',-999)
    base_dir = None  # define prior to running
    # size requirements: one run (one well) is ~ 47 MB there are one full grid run is ~85 GB
    # run time requirments:
    fluxes = [-5, -25, -100]
    fluxes = [e * 86.4 for e in fluxes]  #  define 3 fluxes to run
    model_id = 'opt'  # old re-define
    for flux in fluxes:
        path = os.path.join(base_dir, 'flux_{}'.format(flux))
        well_by_well_depletion_grid(model_id, flux, path, 'one flux of 3, using flux of {} m3/d'.format(flux))
