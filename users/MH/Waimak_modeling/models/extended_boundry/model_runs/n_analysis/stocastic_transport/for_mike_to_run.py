# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/05/2018 4:52 PM
"""

from __future__ import division
from core import env
from gmp_plus_reductions import setup_run_gmp_plus, extract_receptor_data
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

if __name__ == '__main__':
    extract_rep_data = True
    #do not change!!!
    gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")

    # duplicate from here if you need to run a new scenario
    run_mt3d_models = True
    # todo update outpath for netcdf
    # path to the results netcdf
    first_sol_nc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\solutionspackage_.nc")
    if run_mt3d_models:
        print('#### starting gmp plus 25% reduction ####')
        con_shp_path = None  # todo update
        con_att_name = None  # todo update
        rch_con = smt.shape_file_to_model_array(con_shp_path, con_att_name, alltouched=True,
                                                area_statistics=True, fine_spacing=10, resample_method='average')

        # todo update description
        nc_description = 'a 25 percent reduction on all dairy and dairy support in the waimakariri zone'

        # todo update base dir
        # local drive to hold the temporary model files
        base_mt3d_dir = r"D:\mh_waimak_models\base_for_GMP_plus25per_dairy"

        setup_run_gmp_plus(rch_con, base_mt3d_dir, first_sol_nc, nc_description, dt0=1e4, ttsmax=1e5)


    # for new scens duplicate up to here...
    if extract_rep_data:
        extract_receptor_data(scenario_paths={
            'first_sol': first_sol_nc  # todo define scenario name plus add more if needed
        },
            cbc_paths=gmp_cbc,
            outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                    r"ns and results\ex_bd_va\sol_pkgs\all_scens")

        )
