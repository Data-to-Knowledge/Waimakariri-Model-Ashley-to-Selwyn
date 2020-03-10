"""
coding=utf-8
Author: matth
Date Created: 3/07/2017 3:19 PM
"""

from __future__ import division
import fiona
import rasterio.mask
import glob
import os
import numpy as np

if __name__ == '__main__':

    with fiona.open(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Model Grid\model_layering\shp\avonside_chch_extent.shp", "r") as shapefile:
        features = [feature["geometry"] for feature in shapefile]

        rasters = glob.glob(r"C:\Users\MattH\Downloads\layer_thicknesses\*.tif")

        for rast in rasters:
            print(os.path.basename(rast))
            with rasterio.open(rast) as src:
                out_image, out_transform = rasterio.mask.mask(src, features,
                                                                    crop=True)
            out_image[out_image < -999] = np.nan
            print np.nanmean(out_image)
            print np.nanpercentile(out_image,50)
            print ''

