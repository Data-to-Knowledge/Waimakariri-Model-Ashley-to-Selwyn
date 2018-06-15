from __future__ import division
from users.MH.Waimak_modeling.models.extended_boundry.m_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
from scipy.interpolate import rbf
import netCDF4 as nc
import numpy as np
from core import env
import shutil
import os
import pandas as pd
import flopy
from copy import deepcopy
import itertools
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.extract_data import open_path_file_as_df
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.interzone_n import _np_describe
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.mt3d_wrapper import \
    setup_run_mt3d_mp, get_default_mt3d_kwargs, get_sft_stress_period_data, get_ssm_stress_period_data
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import \
    get_gmp_con_layer, get_gmp_plus_con_layer_by_landuse



def timeit_function():
    test = get_gmp_plus_con_layer_by_landuse(False, False, DairySupport=25, DairyFarm=25)

def timeit_function2():
    pass
