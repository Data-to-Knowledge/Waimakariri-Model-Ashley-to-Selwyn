# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 8:37 AM
"""

from __future__ import division
from core import env
import os
import tempfile
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.percentage_reduction_maps import get_pa_reductions
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

def _layer_to_model_array(self, source_layer, attribute, alltouched=False,
                          area_statistics=False, fine_spacing=10,
                          resample_method='average'):  # todo propogate these to upper level
    """
    hidden function to convert a source layer to a rasterized np array
    :param source_layer: from either function above
    :param attribute: attribute to convert.
    :param area_statistics: boolean if true first burn the raster at a finer detail and do some statistics to
                            group up
    :param fine_spacing: int m of spacing to burn the raster to (limits the polygon overlay),
                         only used if area_statistics is True.  fine spacing will be adusted to be a whole factor of
                         self.grid_space
    :param resample_method: key to define resampling options
                            near:
                                nearest neighbour resampling (default, fastest algorithm, worst interpolation
                                quality).
                            bilinear:
                                bilinear resampling.
                            cubic:
                                cubic resampling.
                            cubicspline:
                                cubic spline resampling.
                            lanczos:
                                Lanczos windowed sinc resampling.
                            average:
                                average resampling, computes the average of all non-NODATA contributing pixels.
                                 (GDAL >= 1.10.0)
                            mode:
                                mode resampling, selects the value which appears most often of all the
                                sampled points. (GDAL >= 1.10.0)
                            max:
                                maximum resampling, selects the maximum value from all non-NODATA
                                contributing pixels. (GDAL >= 2.0.0)
                            min:
                                minimum resampling, selects the minimum value from all non-NODATA
                                contributing pixels. (GDAL >= 2.0.0)
                            med:
                                median resampling, selects the median value of all non-NODATA
                                contributing pixels. (GDAL >= 2.0.0)
                            q1:
                                first quartile resampling, selects the first quartile value of all non-NODATA
                                contributing pixels. (GDAL >= 2.0.0)
                            q3:
                                third quartile resampling, selects the third quartile value of all non-NODATA
                                contributing pixels. (GDAL >= 2.0.0)
    :return:  np array of rasterized data
    """
    from osgeo import gdal, osr  # todo see if I can use tempfile for this
    f = tempfile.NamedTemporaryFile()
    temp_file = f.name
    cols = self.cols
    rows = self.rows
    pixelWidth = pixelHeight = self.grid_space  # depending how fine you want your raster

    x_min, y_min = self.ulx, self.uly - self.grid_space * self.rows

    target_ds = gdal.GetDriverByName('GTiff').Create(temp_file, cols, rows,
                                                     1,
                                                     gdal.GDT_Float64)
    target_ds.SetGeoTransform((x_min, pixelWidth, 0, y_min, 0, pixelHeight))
    band = target_ds.GetRasterBand(1)
    NoData_value = -999999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    if alltouched:
        opt = ["ALL_TOUCHED=TRUE", "ATTRIBUTE={}".format(attribute)]
    else:
        opt = ["ATTRIBUTE={}".format(attribute)]
    gdal.RasterizeLayer(target_ds, [1], source_layer, options=opt)

    target_dsSRS = osr.SpatialReference()
    target_dsSRS.ImportFromEPSG(2193)
    target_ds.SetProjection(target_dsSRS.ExportToWkt())
    target_ds = None
    temp_file = '{}/temp.tif'.format(self.temp_file_dir)
    outdata = gdal.Open(temp_file).ReadAsArray()
    if not area_statistics:
        outdata = np.flipud(outdata)
    outdata[np.isclose(outdata, -999999)] = np.nan

    return outdata


if __name__ == '__main__':
    t = _layer_to_model_array(smt,r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\waimakariri_river\waimakariri.shp"
                              'Id',True)
    print('done')