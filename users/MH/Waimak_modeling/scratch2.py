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
from core.classes.hydro import hydro
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

# lat, lon, layer, obs, weigth? i, j
all_targets = pd.read_csv(env.sci(
    "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_inc_error.csv"),
    index_col=0)
all_targets = all_targets.loc[(all_targets.h2o_elv_mean.notnull()) & (all_targets.row.notnull()) &
                              (all_targets.col.notnull()) & (all_targets.layer.notnull())]

# pull out targets for each layer
min_readings = {0: 5,
                1: 2,
                2: 1,
                3: 1,
                4: 1,
                5: 1,
                6: 1,
                7: 1,
                8: 1,
                9: 1}
all_targets.loc[:, 'weight'] = 1 / all_targets.loc[:, 'total_error_m']

outdata = pd.DataFrame()
for layer in range(smt.layers - 1):  # pull out targets for layers 0-9 layer 10 has no targets
    idx = (all_targets.layer == layer) & (all_targets.readings_nondry >= min_readings[layer])
    outdata = pd.concat((outdata, all_targets.loc[idx, ['nztmx', 'nztmy', 'layer', 'h2o_elv_mean',
                                                        'weight', 'row', 'col', 'total_error_m','readings_nondry']]))
no_flow = smt.get_no_flow()
for i in outdata.index:
    layer, row, col = outdata.loc[i, ['layer', 'row', 'col']].astype(int)
    if no_flow[layer, row, col] == 0:  # get rid of non-active wells
        outdata.loc[i, 'layer'] = np.nan
outdata = outdata.dropna(subset=['layer', 'row', 'col'])
outdata.to_csv(r"C:\Users\MattH\Downloads\possible_targets2.csv")