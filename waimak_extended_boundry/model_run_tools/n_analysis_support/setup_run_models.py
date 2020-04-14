# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/04/2018 1:46 PM
"""

from __future__ import division
import env
import os
import pandas as pd
import numpy as np
import shutil
from glob import glob
import logging
import multiprocessing
import psutil
import time
from waimak_extended_boundry.model_run_tools import zipped_modflow_converged, mt3d_converged, setup_run_mt3d_mp, get_default_mt3d_kwargs
from waimak_extended_boundry.model_run_tools.forward_quanity_support.forward_runs import run_forward_runs
from waimak_extended_boundry.model_run_tools.data_extraction.ucn_netcdf import \
    make_ucn_netcd
from waimak_extended_boundry.model_run_tools.data_extraction.cell_budget_netcdf import make_cellbud_netcdf


def start_process():
    """
    function to run at the start of each multiprocess sets the priority lower
    :return:
    """
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def setup_pc5_ftl_repository(model_ids, ftl_dir, base_modelling_dir, increase_eyre=0):
    """

    :param model_ids:
    :param ftl_dir:
    :param base_modelling_dir:
    :param increase_eyre: a value to increase the eyre river flow in cumics
    :return:
    """
    if not os.path.exists(ftl_dir):
        os.makedirs(ftl_dir)

    if not os.path.exists(base_modelling_dir):
        os.makedirs(base_modelling_dir)

    runs = []
    for mid in model_ids:
        runs.append({
            # pc5 (100% efficency via LSR model irrigation, and reduced pumping (25% less to irrigation wells)
            # only in the waimakariri
            'model_id': mid,
            'name': 'pc5_80',
            'base_dir': os.path.join(base_modelling_dir, 'pc5_80'),
            'cc_inputs': None,
            'pc5': True,
            'pc5_well_reduction': True,
            'pc5_to_waimak_only': True,
            'wil_eff': 1,
            'naturalised': False,
            'full_abs': False,
            'pumping_well_scale': 1,
            'org_efficency': 80,
            'write_ftl': True,
            'increase_eyre': increase_eyre
        })

    outputs = run_forward_runs(runs, base_modelling_dir, 'runs to set up a pc5 mt3d_run')
    outputs = pd.DataFrame(outputs, columns=['name', 'success'])
    successful_runs = outputs.loc[outputs.success == 'converged', 'name']
    outputs.to_csv(os.path.join(base_modelling_dir, 'convergence_record.txt'))
    shutil.copyfile(os.path.join(base_modelling_dir, 'convergence_record.txt'),
                    os.path.join(ftl_dir, 'convergence_record.txt'))
    for nm in successful_runs:
        shutil.copyfile(os.path.join(base_modelling_dir, nm, '{}.ftl'.format(nm)),
                        os.path.join(ftl_dir, '{}.ftl'.format(nm)))


def setup_run_mt3d_suite(base_mt3d_dir, ftl_repo, ssm_crch, ssm_stress_period_data, sft_spd, dt0=None, ttsmax=None,
                         scon=None):
    if not os.path.exists(base_mt3d_dir):
        os.makedirs(base_mt3d_dir)
    runs = []
    ftl_paths = glob(os.path.join(ftl_repo, '*.ftl'))
    for ftl_path in ftl_paths:
        name = os.path.basename(ftl_path).replace('.ftl', '')
        default_mt3d_kwargs = get_default_mt3d_kwargs()
        default_mt3d_kwargs.update({'ftl_path': ftl_path,
                                    'mt3d_name': name,
                                    'mt3d_ws': os.path.join(base_mt3d_dir, name),
                                    'ssm_crch': ssm_crch,
                                    'ssm_stress_period_data': ssm_stress_period_data,
                                    'sft_spd': sft_spd,
                                    'safe_mode': False,
                                    'reduce_str_obs': True,
                                    'simplify': False})
        if scon is not None:
            default_mt3d_kwargs['btn_scon'] = scon
        if dt0 is not None:
            default_mt3d_kwargs['dt0'] = dt0
        if ttsmax is not None:
            default_mt3d_kwargs['ttsmax'] = ttsmax

        runs.append(default_mt3d_kwargs)

    t = time.time()
    multiprocessing.log_to_stderr(logging.DEBUG)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map_async(setup_run_mt3d_mp, runs)
    while not results.ready():
        time.sleep(5 * 60)  # sleep 5 min between printing
    pool_outputs = results.get()
    pool.close()  # no more tasks
    pool.join()
    print('completed {} runs in {} min'.format(len(runs), (time.time() - t) / 60))
    with open(os.path.join(base_mt3d_dir, 'status.txt'), 'w') as f:
        f.writelines(['{}\n'.format(e) for e in pool_outputs])


def extract_data(base_mt3d_dir, outfile, description, nname='mednload', units='g/m3'):
    ucn_paths = np.array(glob(os.path.join(base_mt3d_dir, '*', '*.ucn')))
    conv = [mt3d_converged(e.replace('001.UCN', '.list')) for e in ucn_paths]
    ucn_paths = ucn_paths[conv]
    sobs = [e.replace('001.UCN', '.sobs') for e in ucn_paths]
    nsmc_nums = [int(os.path.basename(e).split('_')[0].replace('NsmcReal', '')) for e in ucn_paths]
    make_ucn_netcd(nsmc_nums=nsmc_nums,
                   ucn_paths={nname: ucn_paths},
                   units=units,
                   description=description,
                   nc_path=outfile,
                   zlib=False,
                   ucn_no_value=-1,
                   sobs={nname: sobs})


def extract_cbc_data(base_modflow_dir, description, nc_path, zlib=False):
    sfo_paths = np.array(glob(os.path.join(base_modflow_dir, '*', '*.sfo')))
    conv = np.array([zipped_modflow_converged(e.replace('.sfo', '.nam')) for e in sfo_paths])
    cbc_paths = np.array([e.replace('.sfo', '.cbc') for e in sfo_paths])[conv]
    nsmc_nums = np.array([int(os.path.basename(e).split('_')[0].replace('NsmcReal', '')) for e in sfo_paths])[conv]
    sfo_paths = sfo_paths[conv]
    make_cellbud_netcdf(nsmc_nums, sfo_paths, cbc_paths, description, nc_path, zlib)


pc5_ftl_repo = env.pc5_ftl_repo
if __name__ == '__main__':
    pass
