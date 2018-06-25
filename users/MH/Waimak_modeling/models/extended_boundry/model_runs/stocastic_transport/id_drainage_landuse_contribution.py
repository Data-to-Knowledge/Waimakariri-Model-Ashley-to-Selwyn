# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/06/2018 2:52 PM
"""

from __future__ import division
from core import env
import geopandas as gpd
import os
from gmp_plus_reductions import setup_run_gmp_plus, extract_receptor_data
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

if __name__ == '__main__':
    run_mt3d = True

    base_nload_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                              'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')

    land_types = {'Arable': 2,
                  'DairyFarm': 1,
                  'DairySupport': 1,
                  'ForestTussock': 0,
                  'Horticulture': 2,
                  'Lifestyle': 3,
                  'NotFarm': 0,
                  'OtherinclGolf': 0,
                  'Pigs': 2,
                  'SheepBeefDeer': 2,
                  'SheepBeefHill': 2,
                  'Unknow': 0}

    soil_classes = {'': 0,
                    'WW': 0,
                    'L': 0,
                    'XL': 0,
                    'D': 0,
                    'Pd': 1,
                    'PdL': 2, #todo should this be poorly drained?
                    'F1': 0,
                    'S3': 0,
                    'M': 0,
                    'VL': 0,
                    'S1': 0,
                    'F3': 0,
                    'S2': 0,
                    'S4': 0,
                    'F2': 0,
                    'O': 0
}

    base_data = gpd.read_file(base_nload_path)
    base_data.replace({'luscen_cat': {
        'Forest-Tussock': 'ForestTussock',
        'Other-inclGolf': 'OtherinclGolf',
        'Sheep-Beef-Deer': 'SheepBeefDeer',
        'SheepBeef-Hill': 'SheepBeefHill',
    }})
    base_data.loc[base_data.MGMSoilCod.isnull(), 'MGMSoilCod'] = ''

    base_data.loc[:,'use_soil'] = base_data.MGMSoilCod.replace(soil_classes)
    base_data.loc[:,'use_ltype'] = base_data.luscen_cat.replace(land_types)

    temp_shp_path = os.path.join(smt.temp_file_dir, 'temp_ncon_shp', 'temp_ncon_shp.shp')
    if not os.path.exists(os.path.dirname(temp_shp_path)):
        os.makedirs(os.path.dirname(temp_shp_path))
    base_data.to_file(temp_shp_path)

    soil_class = smt.shape_file_to_model_array(temp_shp_path, attribute='use_soil', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='mode')
    land_class = smt.shape_file_to_model_array(temp_shp_path, attribute='use_ltype', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='mode')

    # todo make raster layer maps and check them




    sft_kwargs = {'eyre': 0, 'waimak': 0, 'ash_gorge': 0, 'cust': 0, 'cust_biwash': 0, 'ash_tribs': 0,
                  'eyre_mar': 0}
    ssm_kwargs = {'wil_race_con': 0, 'upper_n_bnd_flux_con': 0, 'lower_n_bnd_flux_con': 0,
                  'well_stream_con': 0, 'llrzf_con': 0, 'ulrzf_con': 0, 's_race_con': 0, 'chb_con': 0}

    scens = ['dairy_lowdrain',
             'sbda_lowdrain',
             'dairy_highdrain',
             'sbda_highdrain',
             'lifestyle']
    scen_paths = {}
    for scen in scens:
        scen_paths[scen] = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\{}_fraction.nc".format(scen))

    if run_mt3d:
        # setup 5 runs to ID the differnt components
        # use dairy low drainage
        scen = 'dairy_lowdrain'
        out_nc = scen_paths[scen]
        print('#### staring {}'.format(scen))
        rch_con = None  # todo
        nc_description = 'all dairy on poorly drained soils (denoted as pd soils) were set to a concentration of 1 everything else was set to zero'
        base_mt3d_dir = r"D:\mh_waimak_models\{}".format(scen)
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5,
                           ssm_kwargs=ssm_kwargs, sft_kwargs=sft_kwargs)

        # dairy high drainage
        scen = 'dairy_highdrain'
        out_nc = scen_paths[scen]
        print('#### staring {}'.format(scen))
        rch_con = None  # todo
        nc_description = 'all dairy on well drained soils (denoted as not pd soils) were set to a concentration of 1 everything else was set to zero'
        base_mt3d_dir = r"D:\mh_waimak_models\{}".format(scen)
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5,
                           ssm_kwargs=ssm_kwargs, sft_kwargs=sft_kwargs)

        # sbda high drainage
        scen = 'sbda_highdrain'
        out_nc = scen_paths[scen]
        print('#### staring {}'.format(scen))
        rch_con = None  # todo
        nc_description = (
        'all sheep beef deer arable, pig, horticolteral on well drained soils (denoted as not pd soils) were set to a concentration of 1 everything else was set to zero')
        base_mt3d_dir = r"D:\mh_waimak_models\{}".format(scen)
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5,
                           ssm_kwargs=ssm_kwargs, sft_kwargs=sft_kwargs)

        # sbda low drainage
        scen = 'sbda_lowdrain'
        out_nc = scen_paths[scen]
        print('#### staring {}'.format(scen))
        rch_con = None  # todo
        nc_description = 'all sheep beef deer arable, pig on poorly drained soils (denoted as pd soils) were set to a concentration of 1 everything else was set to zero'
        base_mt3d_dir = r"D:\mh_waimak_models\{}".format(scen)
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5,
                           ssm_kwargs=ssm_kwargs, sft_kwargs=sft_kwargs)

        # lifestyle
        scen = 'lifestyle'
        out_nc = scen_paths[scen]
        print('#### staring {}'.format(scen))
        rch_con = None  # todo
        nc_description = 'all sheep beef deer arable, pig on poorly drained soils (denoted as pd soils) were set to a concentration of 1 everything else was set to zero'
        base_mt3d_dir = r"D:\mh_waimak_models\{}".format(scen)
        setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, dt0=1e4, ttsmax=1e5,
                           ssm_kwargs=ssm_kwargs, sft_kwargs=sft_kwargs)

    # extract all data
    extract_data = True
    if extract_data:
        gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")
        extract_receptor_data(scenario_paths=scen_paths,
                              cbc_paths=gmp_cbc,
                              outdir=(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                      r"ns and results\ex_bd_va\zc_n_sol"
                                      r"s\landuse_drainage_fractions"))
