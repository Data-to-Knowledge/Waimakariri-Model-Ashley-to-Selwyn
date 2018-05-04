# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/05/2018 12:41 PM
"""

from __future__ import division
from core import env
from gmp_plus_reductions import setup_run_gmp_plus, extract_receptor_data
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import \
    get_gmp_con_layer

if __name__ == '__main__':
    run_long_term_reductions = True
    extract_data = True
    data_to_report_format = False
    scenarios = { #todo fill out


    }
    raise NotImplementedError
    if run_long_term_reductions:
        for scen, red_layer_path in scenarios.items():
            out_nc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\long_red_{}_ucn.nc".format(scen))
            print('#### starting long reduction {} ####'.format(scen))
            reduction_layer = None #todo!!!
            rch_con = get_gmp_con_layer()
            rch_con *= reduction_layer

            nc_description = 'reductions applied to meet long scenario {}'.format(scen)
            base_mt3d_dir = r"D:\mh_waimak_models\long_red_{}".format(scen)
            setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5)

    gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")
    if extract_data:
        scen_ncs = {}
        for scen in scenarios.keys():
            scen_ncs[scen] = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\long_red_{}_ucn.nc".format(scen))

        extract_receptor_data(scenario_paths=scen_ncs,
                              cbc_paths=gmp_cbc,
                              outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                      r"ns and results\ex_bd_va\zc_n_sol"
                                      r"s\long_scens"))

    if data_to_report_format: #todo make it easy to pull out the data in the report formats
        raise NotImplementedError