"""
Author: matth
Date Created: 27/04/2017 8:37 AM
"""

from __future__ import division

import timeit

import flopy
import numpy as np
import matplotlib.pyplot as plt
import users.MH.Waimak_modeling.model_tools as mt
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

# lat, lon, layer, obs, weigth? i, j
flopy.mbase.run_model(
    exe_name="P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\models_exes\mt3d_usgs_brioch_comp\mt3d-usgs-1.0.exe",
    namefile="NsmcReal000018_pc5_80.nam", model_ws=r"D:\mh_waimak_models\base_for_pc580_mt3d\NsmcReal000018_pc5_80")
