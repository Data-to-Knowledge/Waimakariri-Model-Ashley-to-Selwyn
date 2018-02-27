# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 24/07/2017 8:19 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from osgeo import gdal


if __name__ == '__main__':
    # convert the rasters to arrays
    l0 = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\layer0.tif")
    l1 = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\v2layer1sc1.tif")
    l2_5 = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\v2layer2to5.tif")
    l6_10 = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\v2layer6to10.tif")


    for layer, name in zip([l0,l1,l2_5,l6_10], ['0','1','2-5','6-10']):
        temp = gdal.Open(layer)
        temp2 = temp.GetRasterBand(1).ReadAsArray()
        assert temp2.shape == (smt.rows,smt.cols)
        out_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_{}.txt".format(name))
        np.savetxt(out_path, temp2)

    # make an array of the convined zone
    temp = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/confined_zone_full_extent.shp".format(smt.sdp),
                                         'Id',True)
    temp[np.isfinite(temp)] = 1
    temp[np.isnan(temp)] = 0
    temp = temp.astype(np.int8)
    out_path = env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\confined_bool_array.txt")
    np.savetxt(out_path,temp)