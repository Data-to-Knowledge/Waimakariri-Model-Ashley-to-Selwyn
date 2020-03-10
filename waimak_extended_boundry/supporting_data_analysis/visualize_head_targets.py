# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 23/07/2017 3:41 PM
"""

from __future__ import division
import env
import pandas as pd
from waimak_extended_boundry import smt


if __name__ == '__main__':
    targets2008 = pd.read_csv(env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_inc_error.csv"))
    targets2008 = targets2008.loc[(targets2008.h2o_elv_mean.notnull()) & (targets2008.row.notnull()) & (targets2008.col.notnull())]
    plot_targets = smt.df_to_array(targets2008,'readings',True)
    for layer in range(smt.layers):
        smt.plt_matrix(plot_targets[layer],no_flow_layer=layer,title='layer: {}'.format(layer),vmax=20)

    print('done')