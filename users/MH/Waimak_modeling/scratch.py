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


    [u'ashley_sh1', u'cust_skewbridge', u'cam_bramleys_s', u'cam_marshes_s',
     u'courtenay_kaiapoi_s', u'kaiapoi_harpers_s', u'kaiapoi_island_s',
     u'northbrook_marshes_s', u'ohoka_island_s', u'saltwater_factory_s',
     u'southbrook_marshes_s', u'taranaki_gressons_s', u'taranaki_preeces_s',
     u'waikuku_sh1_s', u'ash_ash_est', u'ash_est_all', u'cam_end_s',
     u'waikuku_end_s', u'taranaki_end_s', u'saltwater_end_s',
     u'kaiapoi_end']

    [u'cam_youngs',
     u'courtenay_neeves',
     u'greigs_greigs',
     u'northbrook_marsh',
     u'ohoka_island',
     u'saltwater_toppings',
     u'silverstream_neeves',
     u'southbrook_marsh',
     u'taranaki_preeces',
     u'waikuku_waikuku-beach-rd',
     u'ashley_sh1',
     u'cust_oxford',
     u'cust_threlkelds',
     u'n7drain_hicklands',
     u'kaiapoi_nroad',
     u'kaiapoi_mainline']

    flow_data = pd.read_csv(r"K:\mh_modeling\stocastic_forward\condensed_data\raw_data\silverstream_neeves_raw_data.csv",index_col=0)
    wai_data = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_str_data.csv", index_col=0)

    data = pd.merge(flow_data, pd.DataFrame(wai_data.loc[:,u'kaiapoi_island_s']), right_index=True,left_index=True)

    print('done')