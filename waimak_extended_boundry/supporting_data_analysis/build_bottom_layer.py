"""
coding=utf-8
Author: matth
Date Created: 2/07/2017 8:35 AM
"""

from __future__ import division
import geopandas as gpd

if __name__ == '__main__':
    per_to_add = 0.10
    min_to_add = 15
    data2 = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Model Grid\model_layering\shp\grid_for_layer1_bot_keep_safe.shp")
    data2.loc[data2.zone=='man-confined','confin'] = -9999
    data = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Model Grid\model_layering\shp\grid_for_layer1_bot_keep_safe_qa.shp")
    data['org_confined'] = data2.loc[:,'confin']
    data['interp'] = -99999
    for i in data.index:
        zone = data.loc[i,'zone']
        if zone =='waimak':
            val = data.loc[i,'waimak']
            adder = val*per_to_add
            if adder < min_to_add:
                adder = min_to_add
            data.loc[i,'interp'] = val - adder
            data.loc[i,'adder'] = adder
        elif zone == 'selwyn':
            val = data.loc[i,'selwyn']
            adder = val* per_to_add
            if adder < min_to_add:
                adder = min_to_add
            data.loc[i,'interp'] = val - adder
            data.loc[i,'adder'] = adder
        elif 'conf' in zone:
            val = data.loc[i,'confin']
            data.loc[i,'interp'] = val
        elif zone =='man-inland':
            val = data.loc[i,'manual']
            adder = val*per_to_add
            if adder < min_to_add:
                adder = min_to_add
            data.loc[i,'interp'] = val - adder
            data.loc[i,'adder'] = adder
        else:
            raise ValueError('unexpedted zone {}'.format(zone))
    if data.interp.min() < -999:
        raise ValueError('null_values')

    data.to_file(driver = 'ESRI Shapefile', filename= r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Model Grid\model_layering\shp\interp_for_layer1_bot.shp")