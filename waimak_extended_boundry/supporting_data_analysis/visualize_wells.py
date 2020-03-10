"""
coding=utf-8
Author: matth
Date Created: 2/07/2017 4:46 PM
"""

from __future__ import division
import env
import pandas as pd
from waimak_extended_boundry import _mt, _elvdb_calc
from osgeo.gdal import Open as gdalOpen
import numpy as np
from users.MH.Waimak_modeling.supporting_data_path import sdp

def well_depth_base_layer_1(targets):
    top = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/tops.tif".format(sdp)).ReadAsArray()
    top[np.isclose(top, -3.40282306074e+038)] = 0

    bottom_layer1 = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/layering/base_layer_1.tif".format(sdp)).ReadAsArray()
    idx = np.isclose(bottom_layer1, -3.40282306074e+038)
    bottom_layer1[idx] = top[idx] - 10 #set nan values to 10 m thick.  all these are out of the no-flow bound
    elv = _elvdb_calc()
    targets.loc[:,'depth_from_l1'] = np.nan
    targets.loc[:,'depth_l8t'] = np.nan
    for well in targets.index:
        try:
            layer, row, col = _mt.convert_coords_to_matix(targets.loc[well,'x'],targets.loc[well,'y'],targets.loc[well,'z'],elv)
            targets.loc[well, 'row'], targets.loc[well, 'col'] = row, col
            targets.loc[well,'layer'] = layer
            targets.loc[well,'layer1_elv'] = bottom_layer1[row, col]
            targets.loc[well,'layer8_top'] = elv[7, row, col]
            mx, my, mz = _mt.convert_matrix_to_coords(row,col,layer,elv)
            targets.loc[well, 'model_x'] = mx
            targets.loc[well, 'model_y'] = my
            targets.loc[well, 'model_z'] = mz
        except AssertionError:
            targets.loc[well, 'row'], targets.loc[well, 'col'] = np.nan, np.nan
            targets.loc[well, 'layer1_elv'] = np.nan

    targets.loc[:,'depth_from_l1'] = targets.loc[:, 'layer1_elv'] - targets.loc[:,'z']
    targets.loc[:,'depth_l8t'] = targets.loc[:, 'layer8_top'] - targets.loc[:,'z']
    return targets



if __name__ == '__main__':
    targets = pd.read_csv(env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/water_levels_for_wells_1_plus_readings.csv"))
    targets = targets.loc[:,['nztmx', 'nztmy', 'mid_screen_elv','readings']]
    targets = targets.rename(columns={'nztmx': 'x', 'nztmy': 'y', 'mid_screen_elv': 'z'})
    targets.index.names = ['well']
    targets['group'] = 'targets'

    targets = well_depth_base_layer_1(targets)
    targets.to_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\wells_1_plus_from_layer1_v2.csv"))
    print 'done'