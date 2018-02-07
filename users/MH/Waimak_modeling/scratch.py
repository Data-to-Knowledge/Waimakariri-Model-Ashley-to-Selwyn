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

if __name__ == '__main__':

    outdir = r"C:\Users\MattH\Downloads\data_for_patrick"
    #os.makedirs(os.path.join(outdir,'elv_db'))
    #os.makedirs(os.path.join(outdir,'hk'))

    # elv_db
    elv_db = smt.calc_elv_db()
    thickness = elv_db[0:-1] - elv_db[1:]
    for l in range(1,smt.layers):
        smt.array_to_raster(os.path.join(outdir,'elv_db','thickness_layer_{:02d}.tif'.format(l)), thickness[l-1])

    # hk
    m = get_model('NsmcBase')
    hk = m.upw.hk.array
    hk[np.isclose(hk, 1e-15)] = np.nan
    for l in range(smt.layers):
        smt.array_to_raster(os.path.join(outdir,'hk','hk_layer_{:02d}'.format(l+1)), hk[l])
