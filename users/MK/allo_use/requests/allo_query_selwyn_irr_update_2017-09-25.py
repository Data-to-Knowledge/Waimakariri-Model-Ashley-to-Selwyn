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
from core.allo_use.plot import allo_plt, allo_multi_plot
from core.misc import printf

#################################
### Parameters

## import parameters

## query parameters
grp_by = ['date', 'take_type']
cwms_zone = ['Selwyn - Waihora']
swaz = 'all'
allo_col = ['irr_area_e']
crc2 = 'all'
take_type = 'all'
use_type = ['irrigation']
#gwaz = ['Valetta', 'Mayfield-Hinds']
#shp = r'C:\ecan\local\Projects\requests\hinds_MAR\2017-03-08\poly1.shp'

#waps_csv = 'C:/ecan/local/Projects/requests/helen/set1/set1.csv'

#take_type = ['Take Surface Water']
years = 'all'

debug = False
sd_only = False
agg_yr = True

## output parameters
base_path = 'E:/ecan/local/Projects/requests'
name = 'kirsty_j'
date = '2017-09-25'

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
allo_plt(lw, start='2000', export_path=export_fig_path, export_name=use_name, end='2017')
allo_plt(lw, start='1970', cat=['tot_allo'], export_path=export_fig_path, export_name=allo_name)

lw2 = lw.loc(axis=0)[:, 'Take Groundwater']
lw2 = lw.loc(axis=0)[:, 'Take Surface Water']

allo_plt(lw2, start='2000', export_path=export_fig_path, export_name=use_name, end='2017')
allo_plt(lw2, start='1970', cat=['tot_allo'], export_path=export_fig_path, export_name=allo_name, end='2017')

#allo_multi_plot(lw, agg_level=1)


#############################
### Check oddities

#otop10 = w_query(allo_use2, grp_by=['dates'], allo_col=allo_col, years=[2015], export=False)
#
#otop11 = otop2[otop2.usage_m3.notnull()]
#

lw2 = lw[lw.take_type == 'Take Groundwater']
lw2 = lw[lw.take_type == 'Take Surface Water']

lw2.sort_values('allo').iloc[-20:]
lw2[lw2.date.isin(['2017-06-30', '2016-06-30'])].sort_values('allo').iloc[-10:]


a1 = allo[allo.crc == 'CRC160803']
a1.iloc[0]

a1 = allo[allo.crc == 'CRC167218']
a1.iloc[0]

a1 = allo[allo.crc == 'CRC990675.1']
a1.iloc[0]





