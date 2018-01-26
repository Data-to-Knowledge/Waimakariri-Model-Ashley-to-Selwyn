"""
Author: matth
Date Created: 3/04/2017 1:39 PM
"""

from __future__ import division
from core import env
import numpy as np
import gdal
import ogr
import osr
import os
from users.MH.Waimak_modeling.supporting_data_path import temp_file_dir

# set up a temporary place to hold the files

if not os.path.exists(temp_file_dir):
    os.makedirs(temp_file_dir)


def geodb_to_model_array(path, shape_name, attribute, alltouched=False):
    """
    create model array from a shapefile in a geodatabase
    :param path: path to geodatabase
    :param shape_name: name of the shape file
    :param attribute: name of the attribute to convert to matrix
    :return: np array of the matrix
    """
    driver = ogr.GetDriverByName("OpenFileGDB")
    ds = driver.Open(path, 0)
    source_layer = ds.GetLayer(shape_name)

    outdata = _layer_to_model_array(source_layer, attribute, alltouched=alltouched)
    return outdata


def shape_file_to_model_array(path, attribute, alltouched=False, rows=190, cols=365):
    """
    shape file to vistas array
    :param path: path to shapefile
    :param attribute: attribute name to convert
    :return:  np.array of the data in model format
    """
    # open shape file
    source_ds = ogr.Open(path)
    source_layer = source_ds.GetLayer()

    outdata = _layer_to_model_array(source_layer, attribute, alltouched=alltouched, rows=rows,cols=cols)
    return outdata


def _layer_to_model_array(source_layer, attribute, alltouched=False, cols=365, rows=190):
    """
    hidden function to convert a source layer to a rasterized np array
    :param source_layer: from either function above
    :param attribute: attribute to convert.
    :return:  np array of rasterized data
    """

    pixelWidth = pixelHeight = 200  # depending how fine you want your raster

    x_min, y_min = 1512162.53275, 5215083.5772 - 200*rows # was 5177083.5772

    target_ds = gdal.GetDriverByName('GTiff').Create('{}/temp.tif'.format(temp_file_dir), cols, rows, 1,
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

    outdata = gdal.Open('{}/temp.tif'.format(temp_file_dir)).ReadAsArray()
    outdata = np.flipud(outdata)
    outdata[np.isclose(outdata,-999999)] = np.nan
    return outdata


# clean up and debug the below and move to core # also add all touched
def geodb_to_raster(path, shape_name, attribute, out_path, pixelWidth, pixelHeight=None):
    """
    convert layer to a raster and save  assumes NZTM.
    :param path: name of the geodb file
    :param shape_name: name of the shape file
    :param source_layer: source layer (from open shapefile)
    :param attribute: attribute to put in raster
    :param out_path: path to save geotiff
    :param pixelWidth: width of individual pixels (in m)
    :param pixelHeight: height of pixel or None then default to square pixels from pixelWidth
    :return:
    """
    driver = ogr.GetDriverByName("OpenFileGDB")
    ds = driver.Open(path, 0)
    source_layer = ds.GetLayer(shape_name)
    _layer_to_raster(source_layer, attribute, out_path, pixelWidth, pixelHeight=pixelHeight)


def shape_to_raster(path, attribute, out_path, pixelWidth, pixelHeight=None):
    """
    convert shapefile to a raster and save  assumes NZTM.
    :param path: path to shape file
    :param source_layer: source layer (from open shapefile)
    :param attribute: attribute to put in raster
    :param out_path: path to save geotiff
    :param pixelWidth: width of individual pixels (in m)
    :param pixelHeight: height of pixel or None then default to square pixels from pixelWidth
    :return:
    """
    # open shape file
    source_ds = ogr.Open(path)
    source_layer = source_ds.GetLayer()
    _layer_to_raster(source_layer, attribute, out_path, pixelWidth, pixelHeight=pixelHeight)


def _layer_to_raster(source_layer, attribute, out_path, pixelWidth, pixelHeight=None):
    """
    convert layer to a raster and save  assumes NZTM.
    :param source_layer: source layer (from open shapefile)
    :param attribute: attribute to put in raster
    :param out_path: path to save geotiff
    :param pixelWidth: width of individual pixels (in m)
    :param pixelHeight: height of pixel or None then default to square pixels from pixelWidth
    :return:
    """
    if pixelHeight is None:
        pixelHeight = pixelWidth

    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    cols = int((x_max - x_min) / pixelHeight)
    rows = int((y_max - y_min) / pixelWidth)

    target_ds = gdal.GetDriverByName('GTiff').Create(out_path, cols, rows, 1, gdal.GDT_Float64)
    target_ds.SetGeoTransform((x_min, pixelWidth, 0, y_min, 0, pixelHeight))
    band = target_ds.GetRasterBand(1)
    NoData_value = -99
    band.SetNoDataValue(NoData_value)
    band.FlushCache()

    gdal.RasterizeLayer(target_ds, [1], source_layer, options=["ATTRIBUTE={}".format(attribute)])

    target_dsSRS = osr.SpatialReference()
    target_dsSRS.ImportFromEPSG(2193)
    target_ds.SetProjection(target_dsSRS.ExportToWkt())
    # clean up and add to core


if __name__ == '__main__':
    test = shape_file_to_model_array(r"C:\Users\MattH\Downloads\test_river.shp", 'test', alltouched=True)