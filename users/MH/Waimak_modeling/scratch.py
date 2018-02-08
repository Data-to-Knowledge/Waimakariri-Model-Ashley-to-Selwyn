from __future__ import division
import numpy as np
import pandas as pd
import flopy
import glob
import matplotlib.pyplot as plt
import os
from core.ecan_io import rd_sql, sql_db
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import get_wel_spd
from pykrige.ok import OrdinaryKriging as okrig
import geopandas as gpd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.extract_data import open_path_file_as_df
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import get_all_well_row_col
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.base_modflow_wrapper import get_model
import time
from core.spatial.vector import xy_to_gpd, points_grid_to_poly, spatial_overlays

if __name__ == '__main__':
    print('starting intersection')
    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    n_zone_shp = gpd.read_file(n_load_path)
    source_area_shp_path = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\probable\Fernside.shp"
    zone = gpd.read_file(source_area_shp_path)
    geometry = zone['geometry'].iloc[0]

    t = time.time()

    sindex = n_zone_shp.sindex
    possible_matches_index = list(sindex.intersection(geometry.bounds))
    possible_matches = n_zone_shp.iloc[possible_matches_index]
    print('took {} s to find possible matches'.format(time.time() - t))
    t = time.time()
    interest_area = gpd.overlay(possible_matches, zone,
                                how='intersection')
    # this was really slow... do a r-tree spatial: http://geoffboeing.com/2016/10/r-tree-spatial-index-python/
    print('took {} s to find specific matches'.format(time.time() - t))

    t = time.time()
    test = spatial_overlays(n_zone_shp,zone)
    print('took {} s for mikes'.format(time.time() - t))
    print ('done')
