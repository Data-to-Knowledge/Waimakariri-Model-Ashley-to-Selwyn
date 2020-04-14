# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/11/2017 8:54 AM
"""

from __future__ import division
import env
import netCDF4 as nc
import flopy_mh as flopy
import numpy as np
from waimak_extended_boundry.model_run_tools import \
    _get_kstkpers
from waimak_extended_boundry import smt
from warnings import warn
import sys
import datetime
import pandas as pd
from itertools import izip_longest
import psutil
import gc
import os


# todo check this for problematic references, ect
def make_ucn_netcd(nsmc_nums, ucn_paths, units, description, nc_path, zlib, ucn_no_value=-1, sobs=None):
    """
    creates as netcdf file for all of the ucn data
    :param nsmc_nums: list of nsmc numbers
    :param ucn_paths: dictionary {variable name:[ucn paths]} ucn paths same order as nsmc_nums
    :param units: either string or dictionary variable name: units) to be passed to the netcdf variable attribute\
    :param description: a description (str) to pass to the netcdf attribute
    :param nc_path: the path to save the netcdf files
    :param zlib: boolean, if True will compress data
    :param ucn_no_value: the no data value for the ucn file (will be convereted to np.nan)
    :param sobs dictionary {variable name:[ucn paths]} ucn paths same order as nsmc_nums or None if None don't add SFR obs
    :return:
    """
    if not isinstance(ucn_paths, dict):
        raise ValueError('expected dictionary for unc_paths, see doc')
    if sobs is not None and not isinstance(sobs, dict):
        raise ValueError('expected dicionary for sobs')

    if isinstance(sobs, dict):
        if sobs.keys() != ucn_paths.keys():
            raise ValueError('expected same keys for sobs and ucns')
        if [len(e) for e in sobs.items()] != [len(e) for e in ucn_paths.items()]:
            raise ValueError('sobs and ucns must have the same lenght of realisations')

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

    row = nc_file.createVariable('row', 'i1', ('row',), fill_value=-9)
    row.setncatts({'units': 'none',
                   'long_name': 'model row',
                   'comments': '1 indexed',
                   'missing_value': -9})

    col = nc_file.createVariable('col', 'i1', ('col',), fill_value=-9)
    col.setncatts({'units': 'none',
                   'long_name': 'model column',
                   'comments': '1 indexed',
                   'missing_value': -9})

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
    av_mem = psutil.virtual_memory().total - 4e9
    file_size = smt.get_empty_model_grid(True)[np.newaxis, :, :, :].nbytes
    num_files = int(av_mem // file_size)
    # some checks
    for var, paths in ucn_paths.items():
        print('extracting data for {}'.format(var))
        if isinstance(units, str):
            addu = units
        elif isinstance(units, dict):
            addu = units[var]
        else:
            raise ValueError('units must be either dict or string, see doc')

        temp = nc_file.createVariable(var, 'f4', ('nsmc_num', 'layer', 'row', 'col'), fill_value=np.nan,
                                      zlib=zlib)
        temp.setncatts({'units': addu,
                        'long_name': var,
                        'missing_value': np.nan})

        temp_ucn_file = flopy.utils.UcnFile(paths[0])
        kstpkper = _get_kstkpers(temp_ucn_file, rel_kstpkpers=-1)[0]  # get the last kstpkper
        if len(temp_ucn_file.get_kstpkper()) > 1:
            warn('more than one kstpkper for {}, using the last kstpkper which is {}'.format(var, kstpkper))
        for i, group in enumerate(grouper(num_files, paths)):
            print('starting set {} to {} of {} for {}'.format(i * 3, i * 3 + num_files, len(paths), var))
            num_not_nan = pd.notnull(list(group)).sum()
            outdata = np.zeros((num_not_nan, smt.layers, smt.rows, smt.cols), dtype=np.float32) * np.nan
            for j, path in enumerate(group):
                if j % 100 == 0:
                    print('reading {} of {}'.format(j, j + 100, num_files))
                if path is None:
                    continue
                ucn_file = flopy.utils.UcnFile(path)
                temp_out = ucn_file.get_data(kstpkper=kstpkper).astype(np.float32)
                temp_out[np.isclose(temp_out, ucn_no_value)] = np.nan
                outdata[j] = temp_out
            temp[i * num_files:i * num_files + num_not_nan] = outdata
            gc.collect()

    if sobs is None:
        return nc_file
    av_mem = psutil.virtual_memory().total - 4e9
    file_size = smt.get_empty_model_grid(False).nbytes
    num_files = int(av_mem // file_size)
    for var, paths in sobs.items():
        print('extracting data for {}'.format(var))
        if isinstance(units, str):
            addu = units
        elif isinstance(units, dict):
            addu = units[var]
        else:
            raise ValueError('units must be either dict or string, see doc')

        temp = nc_file.createVariable('sobs_{}'.format(var), 'f4', ('nsmc_num', 'row', 'col'), fill_value=np.nan,
                                      zlib=zlib)
        temp.setncatts({'units': addu,
                        'long_name': 'SFR concentarations for {}'.format(var),
                        'missing_value': np.nan})

        # assume that only one kstpkper
        for i, group in enumerate(grouper(num_files, paths)):
            print('starting set {} to {} of {} for {}'.format(i * 3, i * 3 + num_files, len(paths), var))
            num_not_nan = pd.notnull(list(group)).sum()
            outdata = np.zeros((num_not_nan, smt.rows, smt.cols), dtype=np.float32) * np.nan
            for j, path in enumerate(group):
                if j % 100 == 0:
                    print('reading {} of {}'.format(j, j + 100, num_files))
                if path is None:
                    continue
                outdata[j] = _get_sfr_con_map(path)
            temp[i * num_files:i * num_files + num_not_nan] = outdata
            gc.collect()

    return nc_file


def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def _get_sfr_con_map(path):
    data = pd.read_table(path, delim_whitespace=True)
    data.loc[:, 'SFR-CONCENTRATION'] = pd.to_numeric(data['SFR-CONCENTRATION'], 'coerce')
    data.loc[:, 'row'] = data.loc[:, 'SFR-NODE'] + 1000
    data.loc[:, 'col'] = data.loc[:, 'SFR-NODE'] + 1000

    mappers = pd.read_hdf(os.path.join(env.sdp_required, 'sfr_data.hdf'), 'reach_data')

    mappers.loc[:, 'reachID'] += 1000 # to prevent replace mapping errors
    mappers = mappers.set_index('reachID')
    row_map = mappers.loc[:, 'i'].to_dict()
    col_map = mappers.loc[:, 'j'].to_dict()

    data = data.replace({'row': row_map, 'col': col_map})
    outdata = smt.df_to_array(data, 'SFR-CONCENTRATION')
    return outdata
