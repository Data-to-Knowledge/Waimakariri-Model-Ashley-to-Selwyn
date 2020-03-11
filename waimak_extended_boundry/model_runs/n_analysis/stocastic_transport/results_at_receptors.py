# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/05/2018 3:31 PM
"""

from __future__ import division
import env
from waimak_extended_boundry.model_run_tools.n_analysis_support.gmp_plus_reductions import extract_receptor_data

if __name__ == '__main__':
    # run on Gwater-02
    gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")
    gmp_eyre_cbc = r"C:\mh_waimak_model_data\GMP_eyre_mar_cbc.nc"
    cmp_cbc = r"C:\mh_waimak_model_data\post_filter1_budget.nc"

    cmp_nc = r"C:\mh_waimak_model_data\mednload_ucn.nc"
    gmp_nc = r"C:\mh_waimak_model_data\GMP_mednload_ucn.nc"
    out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus25per_dairy_ucn.nc")
    out_nc_15 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus15per_dairy_ucn.nc")
    out_nc_eyre_mar = r"C:\mh_waimak_model_data\GMP_mednload_eyre_mar_ucn.nc"
    out_nc_interzone_8 = r"C:\mh_waimak_model_data\GMP_mednload_ucn_8kg_ha_interzone.nc"
    out_nc_50_red = r"C:\mh_waimak_model_data\GMP_mednload_ucn_50_reduc_interzone.nc"

    get_recptor_data = False
    if get_recptor_data:
        extract_receptor_data(scenario_paths={'cmp': cmp_nc,
                                              'gmp': gmp_nc,
                                              'gmp_15_red': out_nc_15,
                                              'gmp_25_red': out_nc_25,
                                              'gmp_eyre_mar': out_nc_eyre_mar,
                                              'interzone_8': out_nc_interzone_8,
                                              'interzone_50_red': out_nc_50_red
                                              },
                              cbc_paths={'cmp': cmp_cbc,
                                         'gmp': gmp_cbc,
                                         'gmp_15_red': gmp_cbc,
                                         'gmp_25_red': gmp_cbc,
                                         'gmp_eyre_mar': gmp_eyre_cbc,
                                         'interzone_8': gmp_cbc,
                                         'interzone_50_red': gmp_cbc
                                         },
                              outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                             r"ns and results\ex_bd_va\zc_n_sol" #todo I could change the output folder so that It is clear what is what
                                             r"s\all_scens"))