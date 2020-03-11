"""
coding=utf-8
Author: matth
Date Created: 30/06/2017 12:44 PM
"""

from __future__ import division
import  pandas as pd
import Waimak_modeling_non_extended.model_tools as mt

data = pd.read_csv(r"C:\Users\MattH\Downloads\verticle_targets.csv", index_col='Well No')

elv = mt.calc_elv_db()
for i in data.index:
    layer, row, col = mt.convert_coords_to_matix(*data.loc[i,['NZTM_x', 'NZTM_y','MidScreenRL']],elv_db=elv)
    lon, lat, node_elv = mt.convert_matrix_to_coords(row,col,layer, elv_db=elv)
    data.loc[i,'node_elv'] = node_elv

data['node_diff'] = data.loc[:,'node_elv'] - data.loc[:,'MidScreenRL']
data.to_csv(r"C:\Users\MattH\Downloads\verticle_targets_with_node.csv")