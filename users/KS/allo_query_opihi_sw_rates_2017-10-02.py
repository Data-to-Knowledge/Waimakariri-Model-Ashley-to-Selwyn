# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:22:17 2016

@author: MichaelEK
"""

from os import path, makedirs
from pandas import merge, read_csv, DataFrame
from core.allo_use import allo_query, allo_ts_proc, allo_proc, allo_gis_proc
from datetime import date
#from Tkinter import Tk
#from tkFileDialog import askdirectory
#from warnings import filterwarnings
#
#filterwarnings('ignore')
#Tk().withdraw()

#################################
### Parameters

## import parameters

allo_gis_csv = 'S:/Surface Water/shared/base_data/usage/allo_gis.csv'
allo_csv = 'S:/Surface Water/shared/base_data/usage/allo.csv'
allo_loc_export_path = r'S:\Surface Water\shared\GIS_base\vector\allocations\allo_gis.shp'

## query parameters
swaz = ['Temuka', 'Opihi SH1', 'Te Nga Wai to – Te Ana a Wai', 'Opihi Saleyards', 'Opihi Rockwood', 'South Opuha', 'North Opuha']
from_date = '2000-01-01'
to_date = '2017-06-30'

## output parameters
base_path = r'S:\Surface Water\shared\projects\otop\opihi\allocation'
name = 'opihi_sw_rates'
date1 = str(date.today())

export_path = path.join(base_path, name + '_' + date1 +  '.csv')

if not path.exists(path.split(export_path)[0]):
    makedirs(path.split(export_path)[0])

###################################
### Initial processing

## allocation processing
allo = allo_proc(export_path=allo_csv)

## Determine locations based on a single WAP per consent...and stuff...
allo_gis = allo_gis_proc(allo, export_shp=allo_loc_export_path, export_csv=allo_gis_csv)

## Convert allocation to time series
allo_ts_rates = allo_ts_proc(allo_gis, start=from_date, end=to_date, freq='sw_rates')

#################################
### Select the necessary data and export

allo_gis1 = allo_gis[allo_gis.swaz.isin(swaz)]

allo_ts = merge(allo_ts_rates, allo_gis1[['crc', 'take_type', 'allo_block', 'swaz']], on=['crc', 'take_type', 'allo_block'])

allo_ts.to_csv(export_path, index=False)















