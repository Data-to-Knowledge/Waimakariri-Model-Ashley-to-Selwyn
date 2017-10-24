# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/10/2017 11:07 AM
"""

from __future__ import division
from core import env
from pykrige.ok import OrdinaryKriging
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import pickle
import numpy as np
import os
import netCDF4 as nc
import pandas as pd
import datetime
import sys
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_from_streams import \
    get_samp_points_df
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
import matplotlib.pyplot as plt

depths = [10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 225]


def get_mask():
    """
    define the mask for kriging from no flow array should be k,i,j
    :param recalc:
    :return:
    """
    elv_db = smt.calc_elv_db()
    no_flow = smt.get_no_flow()
    no_flow[no_flow < 0] = 0
    xs, ys = smt.get_model_x_y(False)
    y, z, x = np.meshgrid(ys, depths, xs)
    z = elv_db[0] - z
    mask = np.zeros(smt.model_array_shape)
    zidx = np.repeat(get_zone_array_index('waimak')[np.newaxis, :, :], 11, axis=0)
    mask[~zidx] = np.nan
    mask[~no_flow.astype(bool)] = np.nan
    mask = np.isnan(mask)
    # depth no_flow

    return mask


def krig_stream(inputdata, stream):
    # 3d krigging on x,y, depth x,y resolution of 200 m ? (e.g. model grid)
    data = inputdata.loc[inputdata[stream].notnull()]
    grid_x, grid_y = smt.get_model_x_y(False)
    # this returns a shape of z,y,x
    k3d_temp, ss3d_temp = np.zeros((len(depths),len(grid_y),len(grid_x)))*np.nan, np.zeros((len(depths),len(grid_y),len(grid_x)))*np.nan
    all_mask = get_mask()
    for layer in range(smt.layers-1):
        idx = data.layer == layer
        val = data.loc[idx, stream].values
        x = data.loc[idx, 'mx'].values
        y = data.loc[idx, 'my'].values
        ok2d = OrdinaryKriging(x=x, y=y, z=val)
        mask = all_mask[layer]
        k, ss = ok2d.execute('masked', grid_x, grid_y,
                                           mask=mask[:,:],
                                           backend='vectorized')  # this may cause a memory error if so switch to 'loop'

        k3d_temp[layer] = k.data
        k3d_temp[layer][k.mask] = np.nan
        ss3d_temp[layer] = ss.data
        ss3d_temp[layer][ss.mask] = np.nan


    return k3d_temp, ss3d_temp


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
    outfile.createDimension('layer', smt.layers-1)

    # create variables
    depth = outfile.createVariable('layer', 'f8', ('layer',), fill_value=np.nan)
    depth.setncatts({'units': 'none',
                     'long_name': 'layer',
                     'missing_value': np.nan})
    depth[:] = range(smt.layers-1)

    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'latitude',
                   'missing_value': np.nan})
    lat[:] = y

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'longitude',
                   'missing_value': np.nan})
    lon[:] = x

    for site in sites:
        k3d, ss3d = krig_stream(data, site)
        site_ss3d = outfile.createVariable('var_{}'.format(site), 'f8', ('layer', 'latitude', 'longitude'),
                                           fill_value=np.nan)
        site_ss3d.setncatts({'units': 'None',
                             'long_name': 'variance of interpolated stream depletion from {}'.format(site),
                             'missing_value': np.nan})
        site_ss3d[:] = ss3d
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
        if var in ['longitude', 'latitude', 'depth']:
            continue

        if 'var' in var:
            vmin, vmax = None, None
        elif 'sd' in var:
            vmin, vmax = 1, 100  # todo vmin and vmax?
        else:
            raise ValueError('shouldnt get here')

        varoutdir = os.path.join(outdir, var)
        if not os.path.exists(varoutdir):
            os.makedirs(varoutdir)

        temp = np.array(data.variables[var])
        for layer in range(smt.layers-1):
            fig, ax = smt.plt_matrix(temp[layer], vmin=vmin, vmax=vmax, cmap='RdBu',
                                     title='{} for flux: {}'.format(var, flux), base_map=True)
            fig.savefig(os.path.join(varoutdir,'layer_{:2d}_{}_flux_{}.png'.format(layer, var, flux)))
            plt.close(fig)


def plot_relationship_3_fluxes():  # I really don't know how to visualise this maybe hold off?
    raise NotImplementedError


def krig_plot_sd_grid(data_path, outdir):
    nc_path = os.path.join(outdir, 'interpolated_{}.nc'.format(os.path.basename(data_path).replace('.csv', '')))
    extract_all_stream_krig(data_path, nc_path)
    plot_out_dir = os.path.join(outdir, 'plots_{}'.format(os.path.basename(nc_path).replace('.nc','')))
    plot_all_streams_sd(nc_path, plot_out_dir)


if __name__ == '__main__':
    mask = get_mask()
    print('done')
