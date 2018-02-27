# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:34:32 2017

@author: MichaelEK
"""

from core.classes.hydro import hydro
import pytest
from os import path, getcwd

## Parameters
mtypes1 = ['flow', 'flow_m', 'swl', 'river_wl_cont_qc', 'river_wl_cont_raw', 'river_flow_cont_raw', 'river_flow_cont_qc', 'river_flow_disc_qc']
mtypes2 = ['gwl', 'gwl_m', 'aq_wl_cont_raw', 'aq_wl_cont_qc', 'aq_wl_disc_qc']
mtypes3 = ['usage']
mtypes4 = ['atmos_precip_cont_raw', 'atmos_precip_cont_qc']

sites1 = [70105]
sites2 = ['K38/1081']
sites3 = ['J38/0774']
sites4 = [313710]

from_date = '2017-02-01'
to_date = '2017-04-01'

qual_codes = [10, 18, 20, 30, 50, 11, 21]

min_count=6
resample_code='W'


## Tests
@pytest.mark.parametrize('mtypes', mtypes1)
def test_ecan_get_data_river(mtypes):
    h1 = hydro().get_data(mtypes=mtypes, sites=sites1, from_date=from_date, to_date=to_date, qual_codes=qual_codes)
    h1._base_stats_fun()
    assert (len(h1._base_stats) == 1)


@pytest.mark.parametrize('mtypes', mtypes2)
def test_ecan_get_data_aq(mtypes):
    h1 = hydro().get_data(mtypes=mtypes, sites=sites2, from_date=from_date, to_date=to_date, qual_codes=qual_codes)
    h1._base_stats_fun()
    assert (len(h1._base_stats) == 1)


def test_ecan_get_data_usage():
    h1 = hydro().get_data(mtypes=mtypes3, sites=sites3, to_date=to_date, qual_codes=qual_codes)
    h1._base_stats_fun()
    assert (len(h1._base_stats) == 1)


@pytest.mark.parametrize('mtypes', mtypes4)
def test_ecan_get_data_atmos(mtypes):
    h1 = hydro().get_data(mtypes=mtypes, sites=sites4, from_date=from_date, to_date=to_date, qual_codes=qual_codes, min_count=min_count, resample_code=resample_code)
    h1._base_stats_fun()
    assert (len(h1._base_stats) == 1)


















