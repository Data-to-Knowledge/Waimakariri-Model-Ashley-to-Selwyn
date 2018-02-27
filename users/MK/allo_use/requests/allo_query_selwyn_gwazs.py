# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:22:17 2016

@author: MichaelEK
"""

from os import path, makedirs
from numpy import in1d
from matplotlib import pyplot as plt
from pandas import merge, read_csv, DataFrame, read_hdf, concat
from core.allo_use import allo_query, allo_ts_proc, allo_filter
from core.allo_use.plot import allo_plt, allo_multi_plot
from core.misc import printf
from core.ecan_io import rd_sql
from geopandas import read_file

#################################
### Parameters

## import parameters
allo_gis_csv = 'E:/ecan/shared/base_data/usage/allo_gis.csv'

## query parameters
gwaz = ['Rakaia-Selwyn', 'Selwyn-Waimakariri']
take_type = ['Take Groundwater']

end = '2017-06-30'

## output parameters
base_path = 'E:/ecan/local/Projects/requests'
name = 'fouad'
date = '2017-07-17'

export_path = path.join(base_path, name, date, name + '_allo.csv')
export_path1 = path.join(base_path, name, date, name + '_allo_with_waps.csv')
export_path2 = path.join(base_path, name, date, name + '_allo_no_waps.csv')

#if not path.exists(export_path):
#    makedirs(export_path)

export_fig_path = path.join(base_path, name, date, 'figures')

allo_name = name + '_past_allo.png'
use_name = name + '_allo_use.png'

#################################
### Read in allocation and usage data and merge data

allo_gis = read_csv(allo_gis_csv)

#################################
### Query data

## Select only the consents within the Ashley swaz group
allo2 = allo_gis[allo_gis.gwaz.isin(gwaz) & allo_gis.take_type.isin(take_type)]

## Select the consents active during the specified dates
allo_ts_ann = allo_ts_proc(allo2, end=end, freq='A')
allo_ts_ann.rename(columns={'allo': 'ann_allo_m3'}, inplace=True)

## Extract the allocation info from allo
allo3 = merge(allo_ts_ann, allo2[['crc', 'take_type', 'allo_block', 'wap', 'gwaz']], on=['crc', 'take_type', 'allo_block', 'wap'], how='left')

allo4 = allo3.groupby(['gwaz', 'date']).sum()
allo5 = allo4.unstack(0)

#################################
#### Output results

allo5.to_csv(export_path)
allo2.to_csv(export_path1)

num_cols = ['max_rate', 'max_vol', 'daily_vol', 'cav']
qual_cols = ['gwaz', 'status_details', 'from_date', 'to_date', 'take_type', 'in_gw_allo']

allo7_num = allo2.groupby('crc')[num_cols].sum()
allo7_qual = allo2.groupby('crc')[qual_cols].first()

allo7 = concat([allo7_num, allo7_qual], axis=1)

allo7.to_csv(export_path2)

#################################
#### Plot
#
#allo_ts_ann.rename(columns={'allo': 'tot_allo'}, inplace=True)
#
#df = allo_ts_ann.groupby(['date']).sum()
#
#allo_plt(df, start='2004', cat=['tot_allo'], export_path=export_fig_path, export_name=use_name)
#allo_plt(lw, start='1970', cat=['tot_allo'], export_path=export_fig_path, export_name=allo_name)

#############################
### Check oddities

t1 = allo2[allo2.status_details == 'Issued - Active']
t1.cav.sum()

t2 = allo_filter(allo2)
t3 = t2[t2.status_details == 'Issued - Active']
t3.cav.sum()
