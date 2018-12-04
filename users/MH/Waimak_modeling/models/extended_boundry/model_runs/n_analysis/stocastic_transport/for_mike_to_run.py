# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/05/2018 4:52 PM
"""
from __future__ import division
import sys
sys.path.insert(0, r"C:\Users\michaelek\git\Ecan.Science.Python.Base")
# sys.path.append(r"C:\Users\michaelek\git\Ecan.Science.Python.Base\core")
# sys.path.append(r"C:\Users\michaelek\git\Ecan.Science.Python.Base\flopy_mh")
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
    first_sol_nc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\solutions_package1_2018-11-28.nc")
    if run_mt3d_models:
        print('#### starting gmp plus 25% reduction ####')
        con_shp_path = env.gw_met_data(r"mh_modeling\input_n_soil_conc\input_n_conc_2018-11-30.shp")
        con_att_name = 'ZIPA_conc2'  # todo update
        rch_con = smt.shape_file_to_model_array(con_shp_path, con_att_name, alltouched=True,
                                                area_statistics=True, fine_spacing=10, resample_method='average')

        # todo update description
        nc_description = 'Mike run with new nitrate layer'

        # todo update base dir
        # local drive to hold the temporary model files
        base_mt3d_dir = r"D:\mh_waimak_models\solutions_package1_2018-11-28"

        setup_run_gmp_plus(rch_con, base_mt3d_dir, first_sol_nc, nc_description, dt0=1e4, ttsmax=1e5)


    # for new scens duplicate up to here...
    if extract_rep_data:
        extract_receptor_data(scenario_paths={
            'first_sol': first_sol_nc  # todo define scenario name plus add more if needed
        },
            cbc_paths=gmp_cbc,
            outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\sol_pkgs\all_scens")

        )
