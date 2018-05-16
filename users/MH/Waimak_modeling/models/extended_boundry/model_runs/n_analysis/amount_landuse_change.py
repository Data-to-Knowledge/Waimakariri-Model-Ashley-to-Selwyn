# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/04/2018 6:20 PM
"""

from __future__ import division
from core import env
from glob import glob
import pandas as pd
import itertools
import geopandas as gpd
from percentage_reduction_to_shp import contourftoshpfile
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from osgeo import gdal
import os


def get_percentages(tif_path):
    levels = [0, 3, 7, 20, 100]
    lon, lat = smt.get_model_x_y(False)
    data = gdal.Open(tif_path).ReadAsArray()
    outpath = r"C:\Users\MattH\Downloads\contour_to_shp2.shp"
    contourftoshpfile(data, levels, outpath, lon, lat, extend="both", smooth=False)
    data = gpd.read_file(outpath)
    data.loc[:, 'area'] = data.area
    outdata = []
    for _min in [3, 7, 20]:
        try:
            outdata.append(data.loc[data.loc[:, 'min'] == _min, 'area'].iloc[0] / (1000 * 1000))
        except:
            outdata.append(0)
    return outdata


if __name__ == '__main__':
    scenarios = ['least_pain', 'middle_option', 'most_gain']
    mar_names = ['with_pc5pa00', 'without_pc5pa00']
    # with and without and pc5pa
    mar_pc5pa = [True, False]
    # with mode  = 50th and 95th
    # with and without conservative things
    conservative_shps = ['use_mix']

    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_reductions_use_zones_use_scens"
    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_red_mt3d"
    # mixed
    for ider in ['', '_stream_only', '_wdc_wells_only', '_private_wells_only']:
        outdata = pd.DataFrame(index=pd.MultiIndex.from_product([scenarios, mar_names], names=['option', 'scenario']),
                               columns=['Land areas requiring 3-7% N load reduction',
                                        'Land areas requiring 7-20% N load reduction',
                                        'Land areas requiring >20% N load reduction'])
        for scen, mar_name in itertools.product(scenarios, mar_names):
            tif_path = os.path.join(base_dir, '{}_use_mix_{}{}_reduction.tif'.format(scen, mar_name, ider))
            print(tif_path)
            temp = get_percentages(tif_path)
            outdata.loc[scen, mar_name] = temp
        outdata *= 100
        outdata = outdata.astype(int)
        outdata.to_csv(os.path.join(base_dir, 'areas_mixed{}.csv'.format(ider)))
