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

if __name__ == '__main__':

    pw = pd.read_csv("P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\PrivateWellZones_detailed.csv", index_col=0)
    idx = pw.depth>50
    pw.loc[idx, 'zone_2'] = pw.loc[idx,'Zone_1'] + '_deep'
    pw.loc[~idx, 'zone_2'] = pw.loc[~idx,'Zone_1'] + '_shallow'
    pw.to_csv("P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\PrivateWellZones_detailed.csv")


