"""
 Author: Matt Hanson
 Created: 14/04/2020 3:03 PM
 """

from waimak_extended_boundry.model_run_tools.model_setup.realisation_id import get_stocastic_set
import shutil
import env
import numpy as np
from waimak_extended_boundry.extended_boundry_model_tools import smt
from waimak_extended_boundry.model_run_tools.model_bc_data.LSR_arrays import get_lsrm_base_array
from waimak_extended_boundry.model_run_tools.model_setup.mt3d_wrapper import get_sft_stress_period_data, get_ssm_stress_period_data
from waimak_extended_boundry.model_run_tools.model_bc_data.n_load_layers import get_gmp_con_layer
from waimak_extended_boundry.model_run_tools.metadata_managment.cwms_index import get_zone_array_index
from waimak_extended_boundry.model_run_tools.n_analysis_support.setup_run_models import extract_data, extract_cbc_data,pc5_ftl_repo,setup_pc5_ftl_repository, setup_run_mt3d_suite

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
        temp = smt.shape_file_to_model_array(shp_file_path, 'Id', True)
        rch = get_lsrm_base_array('pc5', None, None, None, 'mean')
        kg_ha = 8  # kg/ha/year
        load = smt.grid_space ** 2 * 0.0001 * kg_ha * 1000  # g/year
        rch_con = load * (1 / (rch * smt.grid_space ** 2 * 365))  # g/yr * yr/m3 = g/m3
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
        temp = smt.shape_file_to_model_array(shp_file_path, 'Id', True)
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

    # chch version
    print('starting 8kg/ha chch')
    run_8_kg_ha_chch = False
    if run_8_kg_ha_chch:
        ssm_crch = smt.get_empty_model_grid()
        kg_ha = 8  # kg/ha/year
        idx = get_zone_array_index(['chch'])
        rch = get_lsrm_base_array('pc5', None, None, None, 'mean')
        load = smt.grid_space ** 2 * 0.0001 * kg_ha * 1000  # g/year
        rch_con = load * (1 / ((rch[
                                    idx].mean() * 1000 * 365) / 1000 * smt.grid_space ** 2))  # g/yr * yr/m3 = g/m3 # assume an average of 150 mm/yr
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
