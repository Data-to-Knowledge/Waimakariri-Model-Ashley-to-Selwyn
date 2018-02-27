# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:17:37 2016

@author: MichaelEK

Template script for naturalisation and misc stats.
"""

from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.misc import printf
from core.ts.sw import stream_nat, flow_stats, malf7d
from core.ts.plot import hydrograph_plot

############################################
#### Parameters

sites = [69618, 69607, 69619, 69615, 69614, 69635, 69602, 69644]

pivot = True

export_path = r'S:\Surface Water\shared\projects\otop\naturalisation\nat_opihi1.csv'
export_allo_data = r'S:\Surface Water\shared\projects\otop\naturalisation\nat_opihi_allo_data.csv'

###########################################
#### Run stream naturalization

## Only return naturalised flows
flow_nat = stream_nat(sites, pivot=pivot, export_path=export_path)

## Return the allo/usage time series in addition to the naturalised flows
flow_nat, allo_use = stream_nat(sites, pivot=pivot, return_data=True, export_path=export_path)
allo_use.to_csv(export_allo_data, index=False)

###########################################
#### Run other stats

flow_nat1 = flow_nat['nat_flow']
nat_stats = flow_stats(flow_nat1)
malf, alfs, mis_alfs, min_mis_alfs, min_date_alf = malf7d(flow_nat1, return_alfs=True)

##########################################
#### Plotting hydrographs

plot_data = flow_nat1[[69614]]['2010-01-01':'2015-01-01'].dropna()
p1 = hydrograph_plot(plot_data, x_period='year', time_format='%Y-%m-%d')

##########################################
#### Testing

#sd1a[sd1a.wap == 'J37/0306']

flow_nat['flow'][69614]






















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



