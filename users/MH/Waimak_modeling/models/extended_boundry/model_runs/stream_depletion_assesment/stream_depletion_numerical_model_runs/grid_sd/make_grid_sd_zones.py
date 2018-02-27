# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 2:04 PM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import geopandas as gpd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


def make_sd_zone(nc_file, layers, outpath, cutoffs):
    """
    make a polygon shapefile for each of the different zones
    :param nc_file: teh nc file of interpolated sd values
    :param layer: which layer to calculate for
    :param cutoff: the pecentage to assign for the cutoff e.g. everthing with 50% or more are in the zone
    :return: gpd
    """
    layers = np.atleast_1d(layers)
    cutoffs = np.atleast_1d(cutoffs)
    data = nc.Dataset(nc_file)
    outdata = []
    out_names = []
    out_layers = []
    out_cutoffs = []
    for var in data.variables.keys():
        if var in ['longitude', 'latitude', 'layer', 'crs', 'sd_waimak_swaz']:
            continue
        for cutoff in cutoffs:
            for layer in layers:
                x, y = smt.get_model_x_y()
                val = np.array(data.variables[var])[layer]
                cs = plt.contourf(x, y, val, levels=[cutoff, 1000000])
                for i in range(len(cs.collections)):
                    ps = cs.collections[i].get_paths()
                    for p in ps:
                        v = p.vertices
                        x = v[:, 0]
                        y = v[:, 1]
                        poly = Polygon([(i[0], i[1]) for i in zip(x, y)])
                        outdata.append(poly)
                        out_names.append(var)
                        out_layers.append(layer)
                        out_cutoffs.append(cutoff)

    out_df = gpd.GeoDataFrame({'geometry': outdata, 'stream': out_names, 'layer': out_layers, 'cutoff': out_cutoffs})
    out_df.to_file(outpath)


if __name__ == '__main__':
    import os
    nc_files = [
        env.gw_met_data(r"mh_modeling\data_to_mark\StrOpt\grid_sd\with_ncarpet_drains\StrOpt_data_2017_10_21\interpolated_StrOpt_sd_grid_data_flux_-432.0.nc"),
        env.gw_met_data(r"mh_modeling\data_to_mark\StrOpt\grid_sd\with_ncarpet_drains\StrOpt_data_2017_10_21\interpolated_StrOpt_sd_grid_data_flux_-2160.0.nc"),
        env.gw_met_data(r"mh_modeling\data_to_mark\StrOpt\grid_sd\with_ncarpet_drains\StrOpt_data_2017_10_21\interpolated_StrOpt_sd_grid_data_flux_-8640.0.nc"),
        env.gw_met_data(r"mh_modeling\data_to_mark\StrOpt\grid_sd\without_ncarpet_drains\StrOpt_data_2017_10_27\interpolated_StrOpt_sd_grid_data_flux_-8640.0.nc")
    ]
    cutoffs = [20, 40, 60, 90]
    layers = range(10)
    for nc_file in nc_files:
        outpath = os.path.join(os.path.dirname(nc_file),
                               'interpolated_zones_flux_{}.shp'.format(os.path.basename(nc_file).split('_')[-1].replace('.nc','')))
        make_sd_zone(nc_file=nc_file,
                     outpath=outpath,
                     cutoffs=cutoffs,
                     layers=layers)
