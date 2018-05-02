# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/05/2018 4:52 PM
"""

from __future__ import division
from core import env
from gmp_plus_reductions import setup_run_gmp_plus, extract_receptor_data
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import \
    get_gmp_plus_con_layer

if __name__ == '__main__':
    run_gmp_plus25 = True
    out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus25per_dairy_ucn.nc")
    if run_gmp_plus25:
        print('#### starting gmp plus 25% reduction ####')
        rch_con = get_gmp_plus_con_layer(DairyFarm=25, DairySupport=25)
        nc_description = 'a 25 percent reduction on all dairy and dairy support in the waimakariri zone'
        base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus25per_dairy"
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc_25, nc_description, dt0=1e4, ttsmax=1e5)

    run_gmp_plus15 = True
    out_nc_15 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus15per_dairy_ucn.nc")
    if run_gmp_plus15:
        print('#### starting gmp plus 15% reduction ####')
        rch_con = get_gmp_plus_con_layer(DairyFarm=15, DairySupport=15)
        nc_description = 'a 15 percent reduction on all dairy and dairy support in the waimakariri zone'
        base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus15per_dairy"
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc_15, nc_description, dt0=1e4, ttsmax=1e5)

    get_recptor_data = False #todo test and then run this
    if get_recptor_data:
        extract_receptor_data(scenario_paths={'gmp_15_red': out_nc_15, #todo add all of the other key scenarios
                                              'gmp_25_red': out_nc_25},
                              cbc_paths=env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc"), #todo if add eyre mar... manage this
                              outdir=env.sci("Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                             "ns and results\ex_bd_va\zone_committee_n_solutions\gmp_15_25_red"))
