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

    missing_wells = [
        'L34/0073',
        'M34/0645',
        'BW24/0370',
        'M35/10410',
        'M35/10410',
        'M35/0335',
        'M35/18871',
        'M35/0335',
        'L35/0893',
        'M35/4238',
        'M35/4238', # duplicate
    ]
    print(set(missing_wells))
    print (('org', len(missing_wells)))
    print (('set', len(set(missing_wells)))) # 3 are duplicates

    all_wells = get_all_well_row_col()
    all_wells.loc[set(missing_wells)].to_csv(r"C:\Users\MattH\Downloads\missing_sd_wells.csv")