# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:22:17 2016

@author: MichaelEK
"""

from os import path, makedirs
from numpy import in1d
from matplotlib import pyplot as plt
from pandas import merge, read_csv, DataFrame
from allo_use_fun import w_query
from allo_use_plot_fun import allo_plt, allo_multi_plot
from misc_fun import printf

#################################
### Parameters

## import parameters
series_csv = 'C:/ecan/shared/base_data/usage/usage_takes_series_sw_up_with_cav.csv'
allo_csv = 'C:/ecan/shared/base_data/usage/takes_results_with_cav_alt.csv'
zones_csv = 'C:/ecan/shared/base_data/usage/WAPs_catch_zones_cwms.csv'

allo_cols = ['crc', 'wap', 'take_type', 'irr_area', 'gw_zone', 'use_type']
zone_cols= ['wap', 'NZTMX', 'NZTMY', 'NIWA_catch', 'catchment', 'g_catch_num', 'g_catch', 'catchment_num', 'swaz_grp', 'swaz', 'swaz_num', 'cwms_zone']

## query parameters
cwms_zone = ['Orari-Temuka-Opihi-Pareora']
#swaz = ['Upper Pareora River']
allo_col = ['ann_allo_m3']

crc_csv = 'C:/ecan/local/Projects/otop/usage/SW_crcs2.csv'

#take_type = ['Take Surface Water']
years = [2015]

## output parameters
base_path = 'C:/ecan/local/Projects/otop/usage'
name = 'new_crc'
date = '2016-11-09'

export_fig_path = path.join(base_path, name, date,'figures')
export_path = path.join(base_path, name, date, name + '_allo_use.csv')

allo_name = name + '_past_allo.png'
use_name = name + '_allo_use.png'

if not path.exists(export_fig_path):
    makedirs(export_fig_path)

#################################
### Read in allocation and usage data and merge data

series = read_csv(series_csv)
allo = read_csv(allo_csv)
sw_zones = read_csv(zones_csv).drop(['FID','take_type'], axis=1)
sw_zones.columns = zone_cols
sw_zones['swaz'] = sw_zones.swaz.str.replace('/', '-')
sw_zones['swaz_grp'] = sw_zones.swaz_grp.str.replace('/', '-')

allo_use1 = merge(series, allo, on=['crc', 'wap'])
allo_use2 = merge(allo_use1, sw_zones, on=['wap'])

### Read in input data to be used in the query

crcs1 = read_csv(crc_csv)
crcs1.columns = ['crc', 'ann_vol']
crcs = crcs1.iloc[:,0].unique().tolist()

#################################
### Query data

lw = w_query(allo_use2, grp_by=['dates'], crc=crcs, allo_col=allo_col, export_path=export_path, years=[2015], export=False, debug=True)

lw[['crc', 'wap', 'ann_allo_m3']].sort_values('crc')


lw2 = lw[['crc', 'ann_allo_m3']].groupby(['crc']).sum()

#index1 = otop1.index.levels[1][~in1d(otop1.index.levels[1], out1)].values
#otop2 = otop1.loc[(slice(None), index1, slice(None)), :]

#################################
### Make directories if necessary

#otop1.index.levels[1]
#
#for i in otop1.index.levels[1]:
#    if not path.exists(path.join(export_path, i)):
#        makedirs(path.join(export_path, i))
#    if not path.exists(path.join(export_path, i, swaz_path)):
#        makedirs(path.join(export_path, i, swaz_path))


################################
### plot the summaries

## Singles
allo_plt(lw, start='2004', export_path=export_fig_path, export_name=use_name)
allo_plt(lw, start='1970', cat=['tot_allo'], export_path=export_fig_path, export_name=allo_name)

#allo_multi_plot(otop2, agg_level=1, export_path=export_path, export_name='_' + name4 + ex_name)


#############################
### Check oddities

#otop10 = w_query(allo_use2, grp_by=['dates'], allo_col=allo_col, years=[2015], export=False)
#
#otop11 = otop2[otop2.usage_m3.notnull()]
#
#PARENT_B1_ALT_ID='CRC951874.1'


allo[['crc', 'wap', 'ind_ann_allo']].sort_values('crc')

allo4 = allo.loc[in1d(allo.crc, crcs), ['crc', 'wap', 'ind_ann_allo']]

allo5 = allo4.groupby('crc').sum().reset_index()
allo6 = merge(crcs1, allo5, on='crc', how='left')

allo6.to_csv(export_path, index=False)

allo[allo.crc == 'CRC050150']




























