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
from osgeo import gdal

if __name__ == '__main__':
    n_load_path = env.sci(
        'Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')


    t2 = smt.shape_file_to_model_array(n_load_path,'nload_cmp',True, area_statistics=True)
    t1 = smt.shape_file_to_model_array(n_load_path,'nload_cmp',True)
    smt.plt_matrix(t1, title='orginal')
    smt.plt_matrix(t2, title='with areastats')
    print('done')



