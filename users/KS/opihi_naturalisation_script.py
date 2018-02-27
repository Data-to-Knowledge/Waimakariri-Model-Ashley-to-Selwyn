# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:17:37 2016

@author: MichaelEK

Template script for naturalisation and misc stats.
"""

from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.misc import printf
from core.ts.sw import flow_stats, malf7d
from opihi_naturalisation_function import stream_nat
from core.ts.plot import hydrograph_plot
from core.classes.hydro import hydro, all_mtypes
from __future__ import division
from numpy import logical_and

############################################
#### Parameters

sites = [69618, 69650, 69607, 69619, 69616, 69615, 69661, 69614, 69635, 69642, 69602, 69601, 69645, 69649, 69643, 69683, 69644, 69621]

wap_exclude = {69635: ['J38/0797'], 69619: ['J37/0240', 'J37/0306'], 69661: ['10913'], 69683: ['J38/0549', 'J38/0193', 'J38/0246'], 69684: ['K38/0842', 'K38/0843', 'K38/0267', 'K38/0841', 'K38/2034', 'K38/1783'], 69614: ['J38/0294', 'J38/0762', 'J38/0814'], 69650: ['J38/0814'], 69607: ['J38/0814', 'I37/0057'], 69618: ['I37/0057'], 69650: ['I37/0057']}  # waps to exclude for each site

#69641: ['K38/0299']

pivot = True

export_path = r'S:\Surface Water\shared\projects\otop\naturalisation\nat_opihi_pre_edits.csv'
export_allo_data = r'S:\Surface Water\shared\projects\otop\naturalisation\nat_opihi_data.csv'

###########################################
#### Run stream naturalization

## Only return naturalised flows
#flow_nat = stream_nat(sites, pivot=pivot, export_path=export_path)

## Return the allo/usage time series in addition to the naturalised flows
flow_nat, allo_use = stream_nat(sites, pivot=pivot, return_data=True, export_path=export_path, wap_exclusion_list=wap_exclude)
allo_use.to_csv(export_allo_data, index=False)

###########################################
flow_nat1 = flow_nat['nat_flow']


###########################################
#Fix errors in the output

#Add usage from missing consent SCY700181 - stock water to selected sites for the whole period of record
flow_nat1.loc[:, (69618, 69650, 69607)] = flow_nat1.loc[:, (69618, 69650, 69607)]+0.017

#Consent CRC992124A - irrigation
CRC992124A_condition = logical_and(flow_nat1.index >= '2002-07-24', flow_nat1.index <='2006-02-15', flow_nat1.index.month.isin([10,11,12,1,2,3,4]))

flow_nat1.loc[CRC992124A_condition, (69618, 69650, 69607)] = flow_nat1.loc[CRC992124A_condition, (69618, 69650, 69607)] + 0.0314

#Consent CRC051460 - irrigation
CRC051460_condition = logical_and(flow_nat1.index >= '2005-02-03', flow_nat1.index.month.isin([10,11,12,1,2,3,4]))

flow_nat1.loc[CRC051460_condition, (69614, 69650, 69607)] = flow_nat1.loc[CRC051460_condition, (69614, 69650, 69607)]+0.0331

#Consent SCY700184 - public water supply
flow_nat1.loc[flow_nat1.index <= '2001-10-01', (69614, 69650, 69607)] = flow_nat1.loc[flow_nat1.index <= '2001-10-01', (69614, 69650, 69607)]+0.011

#Consent CRC970069 - irrigation
CRC970069_condition = logical_and(flow_nat1.index >= '1998-06-04', flow_nat1.index <='2000-01-31', flow_nat1.index.month.isin([10,11,12,1,2,3,4]))

flow_nat1.loc[CRC970069_condition, (69614, 69650, 69607)] = flow_nat1.loc[CRC970069_condition, (69614, 69650, 69607)] + 0.25

#Consent CRC091571 - stockwater
CRC091571_csv = r'S:\Surface Water\shared\projects\otop\naturalisation\Water_use_CRC091571_CRC972200.csv'
CRC091571_data = read_csv(CRC091571_csv, index_col='Date', parse_dates=True)

flow_nat1.loc[flow_nat1.index >= '2003-06-05', 69619] = flow_nat1.loc[flow_nat1.index >= '2003-06-05', 69619]+ CRC091571_data.loc['2003-06-05':'2015-06-30', 'Water use']

#Add usage from missing consent to selected sites for part of the record
"""
flow_nat1.loc[flow_nat1.index > '2015-06-22', (69618, 69650, 69607, 69641)] = flow_nat1.loc[flow_nat1.index > '2009-10-09', (69618, 69650, 69607, 69641)]+0.00025
"""
#Add the missing data for Kakahu at Mulivhills
kakahu_mulvihills = hydro().get_data(mtypes='flow', sites=[69645], qual_codes=[10, 18, 20, 30, 50])
kakahu_mulvihills_data = kakahu_mulvihills.sel_ts(pivot=True)

flow_nat1.loc[flow_nat1.index <= '1999-06-11', (69645)] = kakahu_mulvihills_data.loc[kakahu_mulvihills_data.index <= '1999-06-11', (69645)]

###################################################
#### Run other stats

nat_stats = flow_stats(flow_nat1)
malf, alfs, mis_alfs, min_mis_alfs, min_date_alf = malf7d(flow_nat1, return_alfs=True, max_missing=200)

flow_nat1.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\nat_flows_opihi.csv', index=True)
alfs.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\alfs_opihi.csv', index=True)
nat_stats.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\nat_stats_opihi.csv', index=True)

#####################################################
### Run stats for the unmodifed flow in the Opihi at SH1

uf_sh1_csv = r'S:\Surface Water\shared\projects\otop\naturalisation\uf_sh1.csv'
uf_sh1_flow = read_csv(uf_sh1_csv, index_col='date', parse_dates=True, dayfirst=True)

uf_sh1_stats = flow_stats(uf_sh1_flow)
uf_malf, uf_alfs, uf_mis_alfs, uf_min_mis_alfs, uf_min_date_alf = malf7d(uf_sh1_flow, return_alfs=True)

uf_sh1_stats.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\uf_sh1_stats.csv', index=True)
uf_alfs.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\uf_sh1_alfs_opihi.csv', index=True)

##########################################
#### Plotting hydrographs

plot_data = flow_nat1[[69614]]['2010-01-01':'2015-01-01'].dropna()
p1 = hydrograph_plot(plot_data, x_period='year', time_format='%Y-%m-%d')

##########################################
#### Testing

#sd1a[sd1a.wap == 'J37/0306']

flow_nat['flow'][69614]

#To add extra sd into usage in flow_nat
flow_nat.loc[:, ('sd_rate', [69618, 69650, 69607, 69641])] = flow_nat.loc[:, ('sd_rate', [69618, 69650, 69607, 69641])] + 0.017


#To check if a consent is the source data
sd_hdf='S:/Surface Water/shared/base_data/usage/sd_est_all_mon_vol.h5'
sd = read_hdf(sd_hdf)

sd[sd['crc'] == 'CRC020123']


###############################################\

###Compare to flows calculated in 2015

old_nat_flow_csv = 'S:/Surface Water/backups/MichaelE/Projects/otop/flow/rec_flow_nat_v02.csv'
old_nat_flow = read_csv(old_nat_flow_csv, index_col='date', parse_dates=True)

old_nat_stats = flow_stats(old_nat_flow)
old_malf, old_alfs, old_mis_alfs, old_min_mis_alfs, old_min_date_alf = malf7d(old_nat_flow, return_alfs=True)

old_malf.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\old_malf_opihi.csv', index=True)
old_alfs.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\old_alfs_opihi.csv', index=True)
nat_stats.to_csv(r'S:\Surface Water\shared\projects\otop\naturalisation\old_nat_stats_opihi.csv', index=True)

########################################
### Parameters

sd_path = 'C:/ecan/base_data/usage/sd_est_all_mon_vol.csv'
flow_path = 'C:/ecan/Projects/otop/flow/otop_flow_recorders.CSV'
wap_gauge_path = 'C:/ecan/Projects/otop/GIS/tables/otop_wap_gauge_catch.csv'
gauge_path = 'C:/ecan/Projects/otop/GIS/tables/otop_min_flow_gauges.csv'
gauge_links_path = 'C:/ecan/Projects/otop/GIS/tables/gauge_catch_links.csv'

wap_gauge_cols = ['wap', 'GRIDCODE']
gauge_col = 'SITENUMBER'

export = True

export_rec_flow_path = 'C:/ecan/Projects/otop/flow/rec_flow_nat.csv'
export_gauge_flow_path = 'C:/ecan/Projects/otop/flow/gauge_flow_nat.csv'

########################################
### Naturalize specific gauging sites

flow, gauge, flow_nat, gauge_nat = stream_nat(sd_path, wap_gauge_path, gauge_links_path, flow_path, gauge_path, export_rec_flow_path=export_rec_flow_path, export_gauge_flow_path=export_gauge_flow_path)

#######################################
### Testing section

if ("NIWACatchm" is NULL, NULL, to_int(floor("NIWACatchm" )))
if ( "CatchmentG" is NULL, "niwa_catch", "CatchmentG" )

Con(IsNull("str_raster3"), "Canterbury 8m DEM", "Canterbury 8m DEM" - "str_raster3")

plt1 = regplot(sel1['et'], sel1['irrigation'])
plt2 = regplot(sel1['et'], sel1['other'])

usage_path = 'C:/ecan/base_data/usage/usage_takes_mon_series_all_est_SD_with_cav.csv'

t1 = read_csv(usage_path)

len(t1.wap.unique())

tot_use2.wap.unique()
tot_use3.wap.unique()

sum(sd_et_est1.sd_usage_est.isnull())
sum(sd_est_all_mon_vol.sd_usage_est.isnull())

t1 = wap_gauge[wap_gauge.site.isin([70103])].wap.unique()
t2 = sd1[sd1.wap.isin(t1)]


max_y = 50

site = 69607
site = 70103
site = 70105
site = 69618
site = 69602

t1 = flow[site].dropna()
d1 = flow_nat[site].dropna()

d1.plot(ylim=[0, max_y])
t1.plot(ylim=[0, max_y])



