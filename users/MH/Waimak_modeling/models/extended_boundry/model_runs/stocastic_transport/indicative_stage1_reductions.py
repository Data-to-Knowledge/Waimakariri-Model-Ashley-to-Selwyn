# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/06/2018 2:52 PM
"""

from __future__ import division
from core import env
from gmp_plus_reductions import setup_run_gmp_plus, extract_receptor_data
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import \
    get_gmp_plus_con_layer_by_landuse, get_gmp_plus_con_layer_by_load_limit


if __name__ == '__main__':
    # I ignored pc5pa.
    reductions = [10, 30, 20]

    run_mt3d=True
    if run_mt3d:
        for red in reductions:
            # all land #todo this does not make any sense, perhaps all n load above 5kg/ha?
            out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_on_5kgha_ucn.nc".format(red))
            print('#### starting gmp plus {}% reduction on all land with a gmp nload >= 5kg/ha ####'.format(red))
            rch_con = get_gmp_plus_con_layer_by_load_limit(n_load_limit=5, n_reduction=red) #todo talk to zeb about this
            nc_description = 'a {} percent reduction on  all land with a gmp nload >= 5kg/ha in the waimakariri zone'.format(red)
            base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus_{}per_on_5kgha".format(red)
            setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc_25, nc_description, dt0=1e4, ttsmax=1e5)


            # all land greater than 15 kg/ha
            out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_on_15kgha_ucn.nc".format(red))
            print('#### starting gmp plus {}% reduction on all land with a gmp nload >= 15kg/ha ####'.format(red))
            rch_con = get_gmp_plus_con_layer_by_load_limit(n_load_limit=15, n_reduction=red)
            nc_description = 'a {} percent reduction on  all land with a gmp nload >= 15kg/ha in the waimakariri zone'.format(red)
            base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus{}per_on_15kgha".format(red)
            setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc_25, nc_description, dt0=1e4, ttsmax=1e5)


            # only dairy and dairy support
            out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_dairy_ucn.nc".format(red))
            print('#### starting gmp plus {}% reduction on dairy and dairy support ####'.format(red))
            rch_con = get_gmp_plus_con_layer_by_landuse(DairyFarm=red, DairySupport=red)
            nc_description = 'a {} percent reduction on all dairy and dairy support in the waimakariri zone'.format(red)
            base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus{}per_dairy".format(red)
            setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc_25, nc_description, dt0=1e4, ttsmax=1e5)


    # extract all data
    extract_data = True
    scen_paths = {}
    for red in reductions:
        scen_paths['dairy_red{}'.format(red)] = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_dairy_ucn.nc".format(red))
        scen_paths['kgha_15_red{}'.format(red)] = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_on_15kgha_ucn.nc".format(red))
        scen_paths['kgha_5_red{}'.format(red)] = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus_{}per_on_5kgha_ucn.nc".format(red))

    if extract_data:
        gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")
        extract_receptor_data(scenario_paths=scen_paths,
                              cbc_paths=gmp_cbc,
                              outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                             r"ns and results\ex_bd_va\zc_n_sol"
                                             r"s\stage_1_red"))

