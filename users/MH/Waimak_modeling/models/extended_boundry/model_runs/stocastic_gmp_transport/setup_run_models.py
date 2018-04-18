# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/04/2018 1:46 PM
"""

from __future__ import division
from core import env
import os
import pandas as pd
import shutil
from glob import glob
import logging
import multiprocessing
import psutil
import time
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.min_flows_reliability.forward_runs import \
    run_forward_runs
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.mt3d_wrapper import \
    setup_run_mt3d_mp, get_default_mt3d_kwargs, get_sft_stress_period_data, get_ssm_stress_period_data
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results.ucn_netcdf import \
    make_ucn_netcd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import get_gmp_con_layer

def start_process():
    """
    function to run at the start of each multiprocess sets the priority lower
    :return:
    """
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def setup_pc5_ftl_repository(model_ids, ftl_dir, base_modelling_dir):
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
            'write_ftl': True
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


def setup_run_mt3d_suite(base_mt3d_dir, ftl_repo, ssm_crch, ssm_stress_period_data, sft_spd, dt0=None, ttsmax=None):
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
    with open(os.path.join(base_mt3d_dir, 'status.txt'), 'w') as f:
        f.writelines(['{}\n'.format(e) for e in pool_outputs])


def extract_data(base_mt3d_dir, outfile):
    ucn_paths = glob(os.path.join(base_mt3d_dir, '*', '*.ucn'))
    sobs = [e.replace('.UCN', '.sobs') for e in ucn_paths]
    nsmc_nums = [int(os.path.basename(e).split('_')[0].replace('NsmcReal', '')) for e in ucn_paths]
    make_ucn_netcd(nsmc_nums=nsmc_nums,
                   ucn_paths={'mednload': ucn_paths},
                   units='g/m3',
                   description='the n concentration for the gmp load on a gmp flow simulation (see pc5_80)',
                   nc_path=outfile,
                   zlib=False,
                   ucn_no_value=-1,
                   sobs={'mednload': sobs})


if __name__ == '__main__':
    ftl_repo = r"K:\mh_modeling\pc580_ftls"
    setup_ftls = False
    if setup_ftls:
        setup_pc5_ftl_repository(get_stocastic_set(),
                                 ftl_repo,
                                 r"D:\mh_waimak_models\base_for_pc580_ftls")
    run_mt3d = True
    if run_mt3d: #todo test
        ssm_crch = get_gmp_con_layer()
        ssm_spd = get_ssm_stress_period_data()
        sft_spd = get_sft_stress_period_data()
        setup_run_mt3d_suite(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d",
                             ftl_repo=ftl_repo,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e3, #todo I'm going to try to run this faster for at least the test
                             ttsmax=1e4) # todo I'm going to try to run this faster for at least the test
