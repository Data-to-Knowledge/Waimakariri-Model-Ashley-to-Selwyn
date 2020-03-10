# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 5/04/2018 2:32 PM
"""

from __future__ import division
import pandas as pd
import os
from glob import glob


def sumaraize(maindir):
    paths = glob(os.path.join(maindir, "*_grouped_landuse.csv"))
    load_overviews = pd.read_csv(os.path.join(maindir,'load_overviews.csv'), index_col=0)
    for path in paths:
        nm = os.path.basename(path).replace('_grouped_landuse.csv','')
        temp = pd.read_csv(path, index_col=0)
        load_overviews.loc[nm, 'total_pa_N_kg'] = temp.loc['sum','new_PA_total_load_prorata']

    load_overviews.to_csv(os.path.join(maindir,'load_overviews_with_paN.csv'))

if __name__ == '__main__':
    sumaraize(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\pa_rules_second_trance")