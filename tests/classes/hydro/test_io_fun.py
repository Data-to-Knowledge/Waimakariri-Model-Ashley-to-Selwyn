# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:34:32 2017

@author: MichaelEK
"""

from core.classes.hydro import hydro
import pytest
from os import path, getcwd
from geopandas import read_file

py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))
csv_files = ['test_long1.csv', 'test_wide1.csv', 'test_wide2.csv']
param = {'test_long1.csv': {'dformat': 'long', 'time': 'time', 'mtypes': 'mtype', 'sites': 'site', 'values': 'data'}, 'test_wide1.csv': {'dformat': 'wide', 'multiindex': True}, 'test_wide2.csv': {'time': 'time', 'mtypes': 'river_flow_cont_qc', 'dformat': 'wide'}}
extra_csv = 'test_combine.csv'
geo_shp = 'sites_geo.shp'
export_keys = ['mtypes', 'sites']
netcdf1 = 'test_netcdf1.nc'


@pytest.mark.parametrize('csv', csv_files)
def test_io_csv(csv):
    tparam = param[csv].copy()

    ## Read
    h1 = hydro().rd_csv(path.join(py_dir, csv), **tparam)
    h1._base_stats_fun()
    assert (len(h1._base_stats) > 4)

    ## Write
    dformat = tparam['dformat']
    out_param = {}
    if dformat == 'long':
        out_param.update({'pivot': False})
    else:
        out_param.update({'pivot': True})
    h1.to_csv(path.join(py_dir, csv), **out_param)

    ## Read
    h1 = hydro().rd_csv(path.join(py_dir, csv), **tparam)
    h1._base_stats_fun()
    assert (len(h1._base_stats) > 4)


## Base import
tparam = param[csv_files[0]]
h1 = hydro().rd_csv(path.join(py_dir, csv_files[0]), **tparam)
h1._base_stats_fun()
h1_len = len(h1._base_stats)

## Combine test
h2 = hydro().rd_csv(path.join(py_dir, extra_csv), **tparam)
h2._base_stats_fun()
h2_len = len(h2._base_stats)


def test_combine():
    h3 = h1.combine(h2)
    h3._base_stats_fun()
    assert (len(h3._base_stats) == (h1_len + h2_len))

## Test geo import
geo1 = read_file(path.join(py_dir, geo_shp))[['index', 'geometry']]
geo2 = geo1.set_index('index')


def test_add_geo_loc():
    h1.add_geo_loc(geo2)
    assert (len(h1.geo_loc) == 7)


def test_io_netcdf():
    ## Read
    h4 = hydro().rd_netcdf(path.join(py_dir, netcdf1))
    h4._base_stats_fun()
    assert (len(h4._base_stats) > 4)

    ## Write
    h4.to_netcdf(path.join(py_dir, netcdf1))

    ## Read
    h4 = hydro().rd_netcdf(path.join(py_dir, netcdf1))
    h4._base_stats_fun()
    assert (len(h4._base_stats) > 4)













