# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/04/2018 1:46 PM
"""

from __future__ import division
from core import env
import os
import pandas as pd
import numpy as np
import shutil
from glob import glob
import logging
import multiprocessing
import psutil
import time
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.convergance_check import zipped_modflow_converged, mt3d_converged
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.min_flows_reliability.forward_runs import \
    run_forward_runs
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.mt3d_wrapper import \
    setup_run_mt3d_mp, get_default_mt3d_kwargs, get_sft_stress_period_data, get_ssm_stress_period_data
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results.ucn_netcdf import \
    make_ucn_netcd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import get_gmp_con_layer
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results.cell_budget_netcdf import make_cellbud_netcdf
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.LSR_arrays import get_lsrm_base_array
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index

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
            'increase_eyre':increase_eyre
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


def setup_run_mt3d_suite(base_mt3d_dir, ftl_repo, ssm_crch, ssm_stress_period_data, sft_spd, dt0=None, ttsmax=None, scon=None):
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
    print('completed {} runs in {} min'.format(len(runs),(time.time()-t)/60))
    with open(os.path.join(base_mt3d_dir, 'status.txt'), 'w') as f:
        f.writelines(['{}\n'.format(e) for e in pool_outputs])


def extract_data(base_mt3d_dir, outfile, description, nname='mednload', units='g/m3'):
    ucn_paths = np.array(glob(os.path.join(base_mt3d_dir, '*', '*.ucn')))
    conv = [mt3d_converged(e.replace('001.UCN','.list')) for e in ucn_paths]
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
    sfo_paths = np.array(glob(os.path.join(base_modflow_dir,'*','*.sfo')))
    conv = np.array([zipped_modflow_converged(e.replace('.sfo','.nam')) for e in sfo_paths])
    cbc_paths = np.array([e.replace('.sfo','.cbc') for e in sfo_paths])[conv]
    nsmc_nums = np.array([int(os.path.basename(e).split('_')[0].replace('NsmcReal', '')) for e in sfo_paths])[conv]
    sfo_paths = sfo_paths[conv]
    make_cellbud_netcdf(nsmc_nums, sfo_paths, cbc_paths, description, nc_path, zlib)


pc5_ftl_repo = env.gw_met_data("mh_modeling/pc580_ftls")
if __name__ == '__main__':

    # ####PC5/GMP flow ####
    setup_ftls = False  # a method to stop re-run
    if setup_ftls:
        setup_pc5_ftl_repository(get_stocastic_set(),
                                 pc5_ftl_repo,
                                 r"D:\mh_waimak_models\base_for_pc580_ftls")

    extract_cbcs = False
    if extract_cbcs:
        description = 'the cbc for the gmp flow simulation (see pc5_80)'
        extract_cbc_data(base_modflow_dir=r"D:\mh_waimak_models\base_for_pc580_ftls",
                         description=description,
                         nc_path=r"C:\mh_waimak_model_data\GMP_cbc.nc")

    run_mt3d = False
    if run_mt3d:
        # as present it took about 10 hours to run the set
        ssm_crch = get_gmp_con_layer()
        ssm_spd = get_ssm_stress_period_data()
        sft_spd = get_sft_stress_period_data()
        setup_run_mt3d_suite(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d",
                             ftl_repo=pc5_ftl_repo,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e4,  # I'm going to try to run this faster for at least the test
                             ttsmax=1e5)  # I'm going to try to run this faster for at least the test

    condence_gmp_results = False
    if condence_gmp_results:
        description = 'the n concentration for the gmp load on a gmp flow simulation (see pc5_80)'
        extract_data(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d",
                     outfile=r"C:\mh_waimak_model_data\GMP_mednload_ucn.nc",
                     description=description,
                     nname='mednload',
                     units='g/m3')
        shutil.copyfile(r"C:\mh_waimak_model_data\GMP_mednload_ucn.nc",
                        r"K:\mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn.nc")

    # #### gmp flow with 8kg/ha n load in interzone ####
    print('starting 8kg/ha interzone')
    run_8_kg_ha = False
    if run_8_kg_ha:
        shp_file_path = env.sci("Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and resul"
                                "ts\ex_bd_va\capture_zones_particle_tracking\source_zone_polygon"
                                "s\interzone\conservative_interzone.shp")
        ssm_crch = get_gmp_con_layer()
        temp = smt.shape_file_to_model_array(shp_file_path,'Id',True)
        rch = get_lsrm_base_array('pc5',None,None,None,'mean')
        kg_ha = 8  # kg/ha/year
        load = smt.grid_space**2 * 0.0001 * kg_ha *1000  # g/year
        rch_con = load * (1/(rch*smt.grid_space**2*365))  #g/yr * yr/m3 = g/m3
        idx = (np.isfinite(temp)) & (ssm_crch > rch_con)
        ssm_crch[idx] = rch_con[idx]
        ssm_spd = get_ssm_stress_period_data()

        sft_spd = get_sft_stress_period_data()

        setup_run_mt3d_suite(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_8kg_ha",
                             ftl_repo=pc5_ftl_repo,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e4,  # I'm going to try to run this faster for at least the test
                             ttsmax=1e5)  # I'm going to try to run this faster for at least the test

    condence_8kgha_results = False
    if condence_8kgha_results:
        description = ('the n concentration for the gmp load and a 8kg/ha load on interzone on a gmp flow '
                       'simulation (see pc5_80)')
        extract_data(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_8kg_ha",
                     outfile=r"C:\mh_waimak_model_data\GMP_mednload_ucn_8kg_ha_interzone.nc",
                     description=description,
                     nname='mednload',
                     units='g/m3')
        shutil.copyfile(r"C:\mh_waimak_model_data\GMP_mednload_ucn_8kg_ha_interzone.nc",
                        r"K:\mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_8kg_ha_interzone.nc")

    print('starting 50% reduc interzone')
    run_50reduc_ha = True
    if run_50reduc_ha:
        shp_file_path = env.sci("Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and resul"
                                "ts\ex_bd_va\capture_zones_particle_tracking\source_zone_polygon"
                                "s\interzone\conservative_interzone.shp")
        ssm_crch = get_gmp_con_layer()
        temp = smt.shape_file_to_model_array(shp_file_path,'Id',True)
        idx = np.isfinite(temp)
        ssm_crch[idx] *= 0.5
        ssm_spd = get_ssm_stress_period_data()

        sft_spd = get_sft_stress_period_data()

        setup_run_mt3d_suite(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_interzone_50red",
                             ftl_repo=pc5_ftl_repo,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e4,  # I'm going to try to run this faster for at least the test
                             ttsmax=1e5)  # I'm going to try to run this faster for at least the test

    condence_50_reduc_results = True
    if condence_50_reduc_results:
        description = ('the n concentration for the gmp load and a 50% reduction of load on interzone on a gmp flow '
                       'simulation (see pc5_80)')
        extract_data(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_interzone_50red",
                     outfile=r"C:\mh_waimak_model_data\GMP_mednload_ucn_50_reduc_interzone.nc",
                     description=description,
                     nname='mednload',
                     units='g/m3')
        shutil.copyfile(r"C:\mh_waimak_model_data\GMP_mednload_ucn_50_reduc_interzone.nc",
                        r"K:\mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_50_reduc_interzone.nc")

    #chch version
    print('starting 8kg/ha chch')
    run_8_kg_ha_chch = False
    if run_8_kg_ha_chch:
        ssm_crch = smt.get_empty_model_grid()
        kg_ha = 8  # kg/ha/year
        idx = get_zone_array_index(['chch'])
        rch = get_lsrm_base_array('pc5',None,None,None,'mean')
        load = smt.grid_space**2 * 0.0001 * kg_ha *1000  # g/year
        rch_con = load * (1/((rch[idx].mean()*1000*365)/1000*smt.grid_space**2))  #g/yr * yr/m3 = g/m3 # assume an average of 150 mm/yr
        ssm_crch[idx] = rch_con
        ssm_spd = get_ssm_stress_period_data()

        sft_spd = get_sft_stress_period_data()

        setup_run_mt3d_suite(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_8kg_ha_chch",
                             ftl_repo=pc5_ftl_repo,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e4,  # I'm going to try to run this faster for at least the test
                             ttsmax=1e5)  # I'm going to try to run this faster for at least the test

    condence_8kgha_chch_results = False
    if condence_8kgha_chch_results:
        description = ('the n concentration for the gmp load and a 8kg/ha load on the christchurch westmelton zone'
                       ' on a gmp flow '
                       'simulation (see pc5_80)')
        extract_data(base_mt3d_dir=r"D:\mh_waimak_models\base_for_pc580_mt3d_8kg_ha_chch",
                     outfile=r"C:\mh_waimak_model_data\GMP_mednload_ucn_8kg_ha_chch.nc",
                     description=description,
                     nname='mednload',
                     units='g/m3')
        shutil.copyfile(r"C:\mh_waimak_model_data\GMP_mednload_ucn_8kg_ha_chch.nc",
                        r"K:\mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_8kg_ha_chch.nc")

    # ####PC5/GMP flow + Eyre Mar ####

    # run an increased eyre river option (Eyre Mar)
    # set additional water to enter at top of segment 4 which is near oxford
    increase_eyre_ftl = r"K:\mh_modeling\pc580_eyre_mar_ftls"
    increase_eyre_cumics = 1
    increase_eyre_mt3d_base = r"D:\mh_waimak_models\base_for_pc580_mt3d_eyre_mar"

    setup_increase_eyre_ftl = False
    if setup_increase_eyre_ftl:
        setup_pc5_ftl_repository(model_ids=get_stocastic_set(),
                                 ftl_dir=increase_eyre_ftl,
                                 base_modelling_dir=r"D:\mh_waimak_models\base_for_pc580_modflow_eyre_mar",
                                 increase_eyre=2)

    run_mt3d_increase_eyre = False
    if run_mt3d_increase_eyre:
        # as present it took about 10 hours to run the set
        ssm_crch = get_gmp_con_layer()
        ssm_spd = get_ssm_stress_period_data()
        sft_spd = get_sft_stress_period_data(eyre_mar=0.1)
        setup_run_mt3d_suite(base_mt3d_dir=increase_eyre_mt3d_base,
                             ftl_repo=increase_eyre_ftl,
                             ssm_crch={0: ssm_crch},
                             ssm_stress_period_data={0: ssm_spd},
                             sft_spd={0: sft_spd},
                             dt0=1e4,  # I'm going to try to run this faster for at least the test
                             ttsmax=1e5)  # I'm going to try to run this faster for at least the test

    condence_eyre_mar_results = False
    if condence_eyre_mar_results:
        description = ('the n concentration for the gmp load on a gmp flow simulation with an additional {} m3/s '
                       'added to the eyre near oxford (see pc5_80)'.format(increase_eyre_cumics))
        extract_data(base_mt3d_dir=increase_eyre_mt3d_base,
                     outfile=r"C:\mh_waimak_model_data\GMP_mednload_eyre_mar_ucn.nc",
                     description=description,
                     nname='mednload',
                     units='g/m3')
        shutil.copyfile(r"C:\mh_waimak_model_data\GMP_mednload_eyre_mar_ucn.nc",
                        r"K:\mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_eyre_mar_ucn.nc")


    extract_cbcs = False
    if extract_cbcs:
        description = 'the cbc for the gmp flow simulation (see pc5_80) with additional 1 m3/s added to eyre'
        extract_cbc_data(base_modflow_dir=r"D:\mh_waimak_models\base_for_pc580_modflow_eyre_mar",
                         description=description,
                         nc_path=r"C:\mh_waimak_model_data\GMP_eyre_mar_cbc.nc")

