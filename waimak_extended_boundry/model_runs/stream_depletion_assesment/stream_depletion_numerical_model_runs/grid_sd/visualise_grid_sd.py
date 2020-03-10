# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/10/2017 11:07 AM
"""

from __future__ import division
from scipy.interpolate import griddata
from waimak_extended_boundry import smt
import numpy as np
import os
import netCDF4 as nc
import pandas as pd
import datetime
import sys
from waimak_extended_boundry import \
    get_samp_points_df
from waimak_extended_boundry import get_zone_array_index
import matplotlib.pyplot as plt

depths = [10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 225]


def get_mask():
    """
    define the mask for kriging from no flow array should be k,i,j
    :param recalc:
    :return:
    """
    no_flow = smt.get_no_flow()
    no_flow[no_flow < 0] = 0
    mask = np.zeros(smt.model_array_shape)
    zidx = np.repeat(get_zone_array_index('waimak')[np.newaxis, :, :], 11, axis=0)
    mask[~zidx] = np.nan
    mask[~no_flow.astype(bool)] = np.nan
    mask = np.isnan(mask)
    # depth no_flow

    return mask


def interplotate_stream(inputdata, stream):
    # 3d krigging on x,y, depth x,y resolution of 200 m ? (e.g. model grid)
    data = inputdata.loc[inputdata[stream].notnull()]
    grid_x, grid_y = smt.get_model_x_y()
    # this returns a shape of z,y,x
    outdata = smt.get_empty_model_grid(True)*np.nan
    all_mask = get_mask()
    for layer in range(smt.layers - 1):
        idx = data.layer == layer
        val = data.loc[idx, stream].values
        x = data.loc[idx, 'mx'].values
        y = data.loc[idx, 'my'].values
        mask = all_mask[layer]
        xs = grid_x[~mask]
        ys = grid_y[~mask]
        temp = griddata(points=(x, y), values=val, xi=(xs, ys), method='cubic')

        outdata[layer][~mask] = temp

    return outdata[0:smt.layers-1]


def extract_all_stream_krig(data_path, outpath):
    """
    extract all streams from a given run and then save as a netcdf
    :param data_path: path to the extracted relative data
    :param outpath: path to save the nc file
    :return:
    """
    samp_points_df = get_samp_points_df()
    sites = list(samp_points_df[samp_points_df.m_type == 'swaz'].index)
    data = pd.read_csv(data_path, skiprows=1)
    flux = data.loc[:, 'flux'].iloc[0]

    outfile = nc.Dataset(outpath, 'w')
    x, y = smt.get_model_x_y(False)
    # create dimensions
    outfile.createDimension('latitude', len(y))
    outfile.createDimension('longitude', len(x))
    outfile.createDimension('layer', smt.layers - 1)

    # create variables
    depth = outfile.createVariable('layer', 'f8', ('layer',), fill_value=np.nan)
    depth.setncatts({'units': 'none',
                     'long_name': 'layer',
                     'missing_value': np.nan})
    depth[:] = range(smt.layers - 1)

    proj = outfile.createVariable('crs', 'i1') #this works really well... #todo add to core
    proj.setncatts({'grid_mapping_name': "transverse_mercator",
                    'scale_factor_at_central_meridian': 0.9996,
                    'longitude_of_central_meridian': 173.0,
                    'latitude_of_projection_origin': 0.0,
                    'false_easting': 1600000,
                    'false_northing': 10000000,
                    })

    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'latitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_y_coordinate'})
    lat[:] = y

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lon.setncatts({'units': 'NZTM',
                   'long_name': 'longitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_x_coordinate'})
    lon[:] = x

    for site in sites:
        k3d = interplotate_stream(data, site)
        site_k3d = outfile.createVariable('sd_{}'.format(site), 'f8', ('layer', 'latitude', 'longitude'),
                                          fill_value=np.nan)
        site_k3d.setncatts({'units': 'percent of pumping',
                            'long_name': 'stream depletion from {}'.format(site),
                            'missing_value': np.nan})
        site_k3d[:] = k3d

    # set global attributes
    outfile.description = (
        'interpolated stream depletion (and variance) at steady state for flux: {} m3/day'.format(flux))
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(data_path, sys.argv[0])
    outfile.flux = flux
    outfile.flux_units = 'm3/day'
    outfile.close()


def plot_all_streams_sd(nc_path, outdir):
    """
    plot all of the k3d and ss3d for each depth extract_all_stream_krig
    :param nc_path: path to the netcdf file created by
    :param outdir: directory to place everything
    :return:
    """
    data = nc.Dataset(nc_path)
    flux = data.flux
    for var in data.variables.keys():
        if var in ['longitude', 'latitude', 'layer', 'crs']:
            continue

        if 'var' in var:
            vmin, vmax = None, None
        elif 'sd' in var:
            vmin, vmax = 1, 100
        else:
            raise ValueError('shouldnt get here')

        varoutdir = os.path.join(outdir, var)
        if not os.path.exists(varoutdir):
            os.makedirs(varoutdir)

        temp = np.array(data.variables[var])
        for layer in range(smt.layers - 1):
            fig, ax = smt.plt_matrix(temp[layer], vmin=vmin, vmax=vmax, cmap='RdBu', no_flow_layer=layer,
                                     title='{} layer {} for flux: {}'.format(var, layer, flux), base_map=True)
            fig.savefig(os.path.join(varoutdir, 'layer_{:2d}_{}_flux_{}.png'.format(layer, var, flux)))
            plt.close(fig)


def krig_plot_sd_grid(data_path, outdir):
    nc_path = os.path.join(outdir, 'interpolated_{}.nc'.format(os.path.basename(data_path).replace('.csv', '')))
    extract_all_stream_krig(data_path, nc_path)
    plot_out_dir = os.path.join(outdir, 'plots_{}'.format(os.path.basename(nc_path).replace('.nc', '')))
    plot_all_streams_sd(nc_path, plot_out_dir)
    return nc_path


if __name__ == '__main__':
    mask = get_mask()
    print('done')
