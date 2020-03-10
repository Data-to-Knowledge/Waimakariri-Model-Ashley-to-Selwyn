# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/09/2017 10:23 AM
"""

from __future__ import division
from waimak_extended_boundry import smt
import geopandas as gpd
import numpy as np

data = gpd.read_file(r"C:\Users\MattH\Downloads\irrigated_area_1109.shp")
data.loc[:,'type_num'] = data.loc[:,'type']
replace_dict = dict(zip(set(data.loc[:,'type']),range(len(set(data.loc[:,'type'])))))
data = data.replace({'type_num':replace_dict})
data.to_file(r"C:\Users\MattH\Downloads\irrigated_area_1109_num_type.shp",driver='ESRI Shapefile')

no_flow = smt.get_no_flow(0)
irrigation = smt.shape_file_to_model_array(r"C:\Users\MattH\Downloads\irrigated_area_1109_num_type.shp",'type_num',True)
irrigation[no_flow==0] = np.nan

for key, value in replace_dict.items():
    print('{}: {} cells {:2.2f} %'.format(key,np.isclose(irrigation,value).sum(),np.isclose(irrigation,value).sum()/np.isfinite(irrigation).sum()*100))
