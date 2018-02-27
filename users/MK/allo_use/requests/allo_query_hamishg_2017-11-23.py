# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:22:17 2016

@author: MichaelEK
"""

from os import path, makedirs
from numpy import in1d
from matplotlib import pyplot as plt
from pandas import merge, read_csv, DataFrame, read_hdf
from core.allo_use import allo_query
from core.allo_use.plot import allo_plt, allo_multi_plot
from core.misc import printf

#################################
### Parameters

## import parameters

## query parameters
grp_by = ['date']
cwms_zone = 'all'
swaz = 'all'
allo_col = ['allo']
crc_csv = r'E:\ecan\local\Projects\requests\hamish_G\2017-11-23\consents2.csv'
crc2 = read_csv(crc_csv)['crc'].unique().tolist()
take_type = 'all'
use_type = 'all'
#gwaz = ['Valetta', 'Mayfield-Hinds']
#shp = r'C:\ecan\local\Projects\requests\hinds_MAR\2017-03-08\poly1.shp'

#take_type = ['Take Surface Water']
years = 'all'

debug = True
sd_only = False
agg_yr = True

## output parameters
base_path = 'E:/ecan/local/Projects/requests'
name = 'hamish_G'
date = '2017-11-23'

export_fig_path = path.join(base_path, name, date, 'figures')
export_path = path.join(base_path, name, date, name + '_allo_use.csv')

allo_name = name + '_past_allo.png'
use_name = name + '_allo_use.png'

if not path.exists(export_fig_path):
    makedirs(export_fig_path)

#################################
### Read in allocation and usage data and merge data


### Read in input data to be used in the query

#crcs = read_csv(waps_csv).RecordNo.unique().tolist()

#################################
### Query data

lw = allo_query(grp_by=grp_by, swaz=swaz, crc=crc2, cwms_zone=cwms_zone, take_type=take_type, use_type=use_type, allo_col=allo_col, years=years, export_path=export_path, sd_only=sd_only, agg_yr=agg_yr, export=True, debug=debug)


#index1 = otop1.index.levels[1][~in1d(otop1.index.levels[1], out1)].values
#otop2 = otop1.loc[(slice(None), index1, slice(None)), :]

#################################
### Daily usage data

use1 = read_hdf(r'E:\ecan\shared\base_data\usage\usage_daily_all.h5').reset_index()
use1.columns = ['wap', 'day', 'usage_day']

### multiple consents on the same month with the same WAP
lw['dup_wap']  = lw.groupby(['wap', 'date'])['allo'].transform('count')

### Combine use with allo/wap
allo_use1 = merge(lw, use1, on=['wap'])
use2 = (allo_use1.loc[:, 'usage_day'] / allo_use1.loc[:, 'dup_wap']).round(2)
allo_use1.loc[:, 'usage_day'] = use2

### Remove unnecessary column
allo_use2 = allo_use1[['crc', 'take_type', 'allo_block', 'wap', 'day', 'usage_day']].drop_duplicates()

### Export
allo_use2.to_csv(path.join(base_path, name, date, 'usage_hamish.csv'), index=False, header=True)






