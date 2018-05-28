# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/05/2018 1:25 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from osgeo.gdal import Open
import numpy as np
import netCDF4 as nc

if __name__ == '__main__':
    # quick script to clean up the interzone shape set the data in the west to NAN
    clipper_path = r'{}\m_ex_bd_inputs\shp\interzone_area_clipper.shp'.format(smt.sdp)
    outpath = r'K:\mh_modeling\interzone_source_zones\stocastic set\individual_netcdfs\cleaned_interzone_area.tif'

    data = nc.Dataset(r'K:\mh_modeling\interzone_source_zones\stocastic set\individual_netcdfs\all_layers_full_tla.nc')
    data = np.array(data.variables['forw_st_number'])

    clipper = smt.shape_file_to_model_array(clipper_path,'Id',True)
    data[np.isfinite(clipper)] = np.nan
    smt.array_to_raster(outpath,data)
