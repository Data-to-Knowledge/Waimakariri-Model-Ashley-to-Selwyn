# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:22:17 2016

@author: MichaelEK
"""

from os import path, makedirs
from numpy import in1d
from matplotlib import pyplot as plt
from pandas import merge, read_csv, DataFrame
from core.allo_use import allo_query
from core.allo_use.plot import allo_plt, allo_multi_plot, allo_restr_plt
from core.misc import printf

#################################
### Parameters

## import parameters

## query parameters
grp_by = ['date', 'catch_name']
cwms_zone = 'all'
#swaz = ['Temuka', 'Opihi Saleyards', 'North Opuha', 'South Opuha', 'Opihi Rockwood', 'Te Nga Wai to – Te Ana a Wai', 'Opihi SH1']
swaz = 'all'
allo_col = ['allo', 'allo_restr']
crc2 = 'all'
take_type = ['Take Surface Water']
use_type = 'all'
#gwaz = ['Valetta', 'Mayfield-Hinds']
#shp = r'C:\ecan\local\Projects\requests\hinds_MAR\2017-03-08\poly1.shp'

#waps_csv = 'C:/ecan/local/Projects/requests/helen/set1/set1.csv'

agg_catch_csv = r'E:\ecan\shared\projects\otop\opihi\allocation\agg_catch.csv'

#take_type = ['Take Surface Water']
years = 'all'

debug = False
sd_only = False
agg_yr = True
rem_over_usage = True

## output parameters
base_path = 'E:/ecan/local/Projects/requests'
name = 'opihi'
date = '2017-11-23'

export_fig_path = path.join(base_path, name, date, 'figures')
export_path = path.join(base_path, name, date, name + '_allo_use.csv')
export_path_all = path.join(base_path, name, date, name + '_allo_use_all.csv')

allo_name = name + '_past_allo.png'
use_name = name + '_allo_use.png'
use_name2 = name + '_allo_use2.png'

if not path.exists(export_fig_path):
    makedirs(export_fig_path)

#################################
### Read in allocation and usage data and merge data

agg_catch = read_csv(agg_catch_csv)
catch_name = agg_catch.catch_name.tolist()

### Read in input data to be used in the query

#crcs = read_csv(waps_csv).RecordNo.unique().tolist()

#################################
### Query data

lw = allo_query(grp_by=grp_by, swaz=swaz, crc=crc2, cwms_zone=cwms_zone, take_type=take_type, use_type=use_type, allo_col=allo_col, years=years, catch_name=catch_name, sd_only=sd_only, agg_yr=agg_yr, debug=debug, rem_over_usage=rem_over_usage)

all1 = allo_query(grp_by=grp_by, swaz=swaz, crc=crc2, cwms_zone=cwms_zone, take_type=take_type, use_type=use_type, allo_col=allo_col, years=years, catch_name=catch_name, export_path=export_path_all, sd_only=sd_only, agg_yr=agg_yr, export=True, debug=True)

lw1 = lw.reset_index()
lw2 = merge(lw1, agg_catch, on='catch_name').drop('catch_name', axis=1)
lw3 = lw2.groupby(['date', 'catchment']).sum()
lw3.to_csv(export_path)
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

allo_multi_plot(lw3, export_path=export_fig_path, export_name=use_name, plot_fun=allo_restr_plt, start='2004', cat=['meter_allo', 'meter_usage'])
allo_multi_plot(lw3, export_path=export_fig_path, export_name=use_name2, plot_fun=allo_restr_plt, start='2004', cat=['tot_allo', 'meter_allo', 'meter_usage'])

allo_restr_plt(lw3, start='2004', cat=['meter_allo', 'meter_usage'], export_path=export_fig_path, export_name=use_name)


#############################
### Check oddities

#otop10 = w_query(allo_use2, grp_by=['dates'], allo_col=allo_col, years=[2015], export=False)
#
#otop11 = otop2[otop2.usage_m3.notnull()]
#













