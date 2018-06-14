# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/11/2017 8:51 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import flopy_mh as flopy
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_at_wells import \
    _get_kstkpers, hds_no_data
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from warnings import warn
import datetime
import sys
import psutil
import pandas as pd
import gc
from itertools import izip_longest

def make_hds_netcdf(nsmc_nums, hds_paths, description, nc_path, zlib):
    """
    make a cell budget file netcdf for easy use
    :param nsmc_nums: list the unique identifiers for the netcdfs
    :param hds_paths: list the paths for the hds files of a given netcdf number (same order as nsmc_nums)
    :param nc_path: the path to the outfile
    :param description: a description (str) to set as a netcdf attribute
    :return:
    """
    print('starting to make a netcdf of heads, this will take a long time')
    nc_file = nc.Dataset(nc_path, 'w')

    # make dimensions
    nc_file.createDimension('nsmc_num', len(nsmc_nums))
    nc_file.createDimension('layer', smt.layers)
    nc_file.createDimension('row', smt.rows)
    nc_file.createDimension('col', smt.cols)

    # general attributes
    nc_file.description = description
    nc_file.history = 'created {}'.format(datetime.datetime.now().isoformat())
    nc_file.source = 'script: {}'.format(sys.argv[0])

    # set up layer col and row
    layer = nc_file.createVariable('layer', 'i1', ('layer',), fill_value=-9)
    layer.setncatts({'units': 'none',
                     'long_name': 'model layer',
                     'comments': '1 indexed',
                     'missing_value': -9})
    layer[:] = range(1, smt.layers + 1)

    row = nc_file.createVariable('row', 'i1', ('row',), fill_value=-9)
    row.setncatts({'units': 'none',
                   'long_name': 'model row',
                   'comments': '1 indexed',
                   'missing_value': -9})
    row[:] = range(1, smt.rows + 1)

    col = nc_file.createVariable('col', 'i1', ('col',), fill_value=-9)
    col.setncatts({'units': 'none',
                   'long_name': 'model column',
                   'comments': '1 indexed',
                   'missing_value': -9})
    col[:] = range(1, smt.cols + 1)

    nsmc_num = nc_file.createVariable('nsmc_num', 'i2', ('nsmc_num',), fill_value=-9)
    nsmc_num.setncatts({'units': 'none',
                        'long_name': 'Null Space Monte Carlo Realisation Number',
                        'comments': 'unique identifier phi lower and phi upper are -1 and -2, respectively',
                        'missing_value': -9})
    nsmc_num[:] = nsmc_nums

    x, y = smt.get_model_x_y(False)

    proj = nc_file.createVariable('crs', 'i1')  # this works really well...
    proj.setncatts({'grid_mapping_name': "transverse_mercator",
                    'scale_factor_at_central_meridian': 0.9996,
                    'longitude_of_central_meridian': 173.0,
                    'latitude_of_projection_origin': 0.0,
                    'false_easting': 1600000,
                    'false_northing': 10000000,
                    })

    lat = nc_file.createVariable('latitude', 'f8', ('row',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'latitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_y_coordinate'})
    lat[:] = y

    lon = nc_file.createVariable('longitude', 'f8', ('col',), fill_value=np.nan)
    lon.setncatts({'units': 'NZTM',
                   'long_name': 'longitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_x_coordinate'})
    lon[:] = x

    av_mem = psutil.virtual_memory().total -4e9
    file_size = smt.get_empty_model_grid(True)[np.newaxis,:,:,:].nbytes
    num_files = int(av_mem//file_size)

    # get the last kstpkper
    temp = flopy.utils.HeadFile(hds_paths[0])
    kstpkper = _get_kstkpers(temp, rel_kstpkpers=-1)[0]
    if len(temp.get_kstpkper()) > 1:
        warn('more than one kstpkper, using the last kstpkper which is {}'.format(kstpkper))

    # create variable
    dim = ('nsmc_num', 'layer', 'row', 'col')
    all_hds = nc_file.createVariable('heads', 'f8', dim, fill_value=np.nan, zlib=zlib)
    all_hds.setncatts({'units': 'm',
                       'long_name': 'modelled heads',
                       'missing_value': np.nan})

    # add data
    for i, group in enumerate(grouper(num_files, hds_paths)):
        print('starting set {} to {} of {}'.format(i * num_files, (i+1) * num_files, len(hds_paths)))
        num_not_nan = pd.notnull(list(group)).sum()
        outdata = np.zeros((num_not_nan, smt.layers, smt.rows, smt.cols), dtype=np.float32) * np.nan
        for j, path in enumerate(group):
            if j % 100 == 0:
                print('reading {} to {} of {}'.format(j, j + 100, num_files))
            if path is None:
                continue
            hds = flopy.utils.HeadFile(path)
            temp_out = hds.get_data(kstpkper=kstpkper).astype(np.float32)
            temp_out[np.isclose(temp_out, hds_no_data)] = np.nan
            outdata[j] = temp_out
        all_hds[i * num_files:i * num_files + num_not_nan] = outdata
        gc.collect()

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)