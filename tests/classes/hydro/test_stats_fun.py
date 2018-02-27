# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:34:32 2017

@author: MichaelEK
"""

from core.classes.hydro import hydro
import pytest
from os import path, getcwd

py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))
netcdf1 = 'test_netcdf2.nc'

## Read in data
h1 = hydro().rd_netcdf(path.join(py_dir, netcdf1))
h1._base_stats_fun()
stats = h1._base_stats


def test_gen_stats():
    flow_stats = h1.stats(mtypes='river_flow_cont', sites=[65101, 69505], below_median=True)
    assert (len(flow_stats) == 2) & (len(flow_stats.columns) == 15)


def test_resample():
    res1 = h1.resample('W')
    res1._base_stats_fun()
    assert (stats['count'][0] > res1._base_stats['count'][0])


def test_malf7d():
    malf, alf, alf_miss, n_days_min, min_date = h1.malf7d(return_alfs=True)
    assert (len(malf) == 5)



