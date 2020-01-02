"""
Author: matth
Date Created: 14/06/2017 3:35 PM
"""

from __future__ import division
from core import env
import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import matplotlib.pylab as plt
from get_str_rch_values import get_base_rch
#from read_write import no_flow

def array_to_raster(path, array, only_active=True):
    if only_active:
        array[no_flow] = -99
    output_raster = gdal.GetDriverByName('GTiff').Create(path, array.shape[1], array.shape[0], 1,
                                                         gdal.GDT_Float32)  # Open the file
    geotransform = (1512162.53275, 200, 0, 5214983.5772000002, 0, -200)
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()  # Establish its coordinate encoding
    srs.ImportFromEPSG(2193)  # This one specifies NZTM.
    output_raster.SetProjection(srs.ExportToWkt())  # Exports the coordinate system
    # to the file
    band = output_raster.GetRasterBand(1)
    band.WriteArray(array)  # Writes my array to the raster
    band.FlushCache()
    band.SetNoDataValue(-99)

if __name__ == '__main__':
    test = get_base_rch()
    array_to_raster(r"C:\Users\MattH\Downloads\test_raster.tif", test)