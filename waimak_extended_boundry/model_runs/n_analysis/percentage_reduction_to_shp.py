# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/04/2018 1:16 PM
"""

from __future__ import division
import os
import numpy as np
from osgeo import ogr, gdal
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import shapely.geometry as shp
from waimak_extended_boundry import smt


def contourftoshpfile(data, levels, out_path, lon, lat, extend="both", smooth=False):
    """
    saves a polygon shape file from a contourf output of the data

    :param data: the array (from chris's output files)
    :param levels: the break points for the dataset
    :param extend: min,max,both,neither.  used to set the extend function of contourf
    :param out_path: path to save the shape file without the file ending (.shp)
    :param lon: extent and number of datapoints for the long and lat grid.  i.e. lon, lat = np.linspace(166.425, 178.525, 729), np.linspace(-34.375,-47.325, 780) #new lat, lon for interpolation
    :param lat: extent and number of datapoints for the long and lat grid.  i.e. lon, lat = np.linspace(166.425, 178.525, 729), np.linspace(-34.375,-47.325, 780) #new lat, lon for interpolation
    :param smooth: bool if true smooth the contours
    """

    out_path = out_path.replace('.shp','')
    # plot contourf
    if smooth:
        # fill NAN space so that there is not blocky shorelines
        mask = np.isnan(data)
        data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
        data = ndimage.interpolation.zoom(np.ma.masked_invalid(data), 3, order=1,
                                                mode='wrap')  # use bilinear interpolation to smooth
        data = np.round(data, 3)
    x, y = np.meshgrid(lon, lat)
    cs = plt.contourf(x, y, data, latlon=True, levels=levels, extend=extend)

    # produce polygons coordinates from contourf
    polygons = list()
    levels2 = [float(l) for l in cs.levels]
    levels2.append(levels2[0]) # add the second one again so it shows up
    for m, polygon in enumerate(cs.collections):
        mpoly = []
        for path in polygon.get_paths():
            path.should_simplify = False
            poly = path.to_polygons()
            exterior = []
            holes = []
            if len(poly) > 0:
                exterior = poly[0]  # and interiors (holes) are in poly[1:]
                # Crazy correction of one vertice polygon, mpl doesn't care about it
                if len(exterior) < 2:
                    continue
                p0 = exterior[0]
                # exterior = np.vstack((exterior, self.epsi_point(p0), self.epsi_point(p0)))#not sure what this is for as I pulled this code from a plugin.
                if len(poly) > 1:  # There's some holes
                    for h in poly[1:]:
                        if len(h) > 2:
                            holes.append(h)

                mpoly.append([exterior, holes])
        if len(mpoly) > 0:
            polygons.append([m, levels2[m-1], levels2[m],
                             mpoly])  # I subtracted 1 becasue the first polygon will never have anything in it effectivly this make the value the min value for the range rather than the max.
        print m, polygon

    # save a shp file
    driver = ogr.GetDriverByName('Esri Shapefile')
    # remove any instances of the same shapefile in the directory so as not to pass on atttributes
    # if the shape file has been written previously this will throw an error
    # to avoid at present I simply restart QGIS
    try:
        os.remove(out_path + ".shp")
    except:
        pass

    try:
        os.remove(out_path + ".dbf")
    except:
        pass

    try:
        os.remove(out_path + ".shx")
    except:
        pass

    ds = driver.CreateDataSource(out_path + ".shp")
    layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('min', ogr.OFTReal))  # add a second attribute for the level
    layer.CreateField(ogr.FieldDefn('max', ogr.OFTReal))  # add a second attribute for the level
    layer.CreateField(ogr.FieldDefn('midpoint', ogr.OFTReal))  # add a second attribute for the level
    defn = layer.GetLayerDefn()
    # add polygon data to shape file
    for l, min, max, polygon in polygons:
        poly = shp.MultiPolygon(polygon)
        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(defn)
        feat.SetField('id', l)
        feat.SetField('min', min)  # working here
        feat.SetField('max', max)  # working here
        feat.SetField('midpoint', (min + max)/2)  # working here

        # Make a geometry, from Shapely object
        try:
            geom = ogr.CreateGeometryFromWkb(poly.wkb)
        except ValueError:
            continue
        feat.SetGeometry(geom)

        layer.CreateFeature(feat)
        feat = geom = None  # destroy these

    # Save and close everything
    ds = layer = feat = geom = None
    return


if __name__ == '__main__':
    levels = [0, 3, 7, 20, 100]
    lon,lat = smt.get_model_x_y(False)
    temp_file = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_reductions_use_zones_excl_interzone_3scen\most_gain_load_95th_use_mix_with_pc5pa00_private_wells_only_reduction.tif"
    outpath = r"C:\Users\MattH\Downloads\testcontour_to_shp2"
    data = gdal.Open(temp_file).ReadAsArray()

    contourftoshpfile(data, levels, outpath, lon, lat, extend="both", smooth=False)