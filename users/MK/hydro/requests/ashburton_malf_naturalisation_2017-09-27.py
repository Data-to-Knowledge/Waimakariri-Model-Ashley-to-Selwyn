# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas for MALf calcs.
"""

#### Import hydro class
from core.classes.hydro import hydro
from os.path import join
from core.ts.sw import stream_nat, flow_stats, malf7d
from core.spatial import rec_catch_del

#################################################
#### Load data

### Parameters

base_path = r'E:\ecan\local\Projects\requests\adam\2017-09-27'

## Import
mtypes1 = 'flow'
sites1 = [68829, 68810, 68819, 68822, 68806, 68801]
qual_codes = [10, 18, 20, 50]

mtypes2 = 'flow_m'
sites2 = [168834, 168821]

bound_shp = 'ashburton.shp'
sites_shp = 'min_flow_sites.shp'

## Export
stats_csv = 'rec_stats.csv'
malf_csv = 'rec_malf.csv'
alf_csv = 'rec_alf.csv'
days_mis_csv = 'rec_alf_days_mis.csv'
reg_csv = 'min_flow_reg.csv'
catch_del_shp = 'min_flow_catch_del.shp'
min_flow_nat_csv = 'min_flow_nat.csv'
rec_flow_nat_csv = 'rec_flow_nat.csv'
min_flow_allo_use_csv = 'min_flow_allo_use.csv'
rec_flow_allo_use_csv = 'rec_flow_allo_use.csv'
min_flow_malf_csv = 'min_flow_malf.csv'
rec_nat_malf_csv = 'rec_nat_malf.csv'
min_flow_nat_malf_csv = 'min_flow_nat_malf.csv'

## Plotting
start = '1986-07-01'
end = '1987-06-30'
x_period = 'month'
time_format = '%d-%m-%Y'

#flow_sites = 70105

################################################
#### Import data

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)
h2 = hydro().get_data(mtypes=mtypes2, sites=sites2)
h3 = h2.get_data(mtypes='flow', sites=join(base_path, bound_shp), qual_codes=qual_codes)

################################################
#### Tools

### Flow tools

## MALF and flow stats
fstats = h1.stats(mtypes=mtypes1)
fstats.to_csv(join(base_path,stats_csv))

#malf1 = h1.malf7d()
#malf1

# Get out ALFs and a whole lot of other things
malf4, alfs, mis_alfs, min_mis_alfs, min_date_alf = h1.malf7d(return_alfs=True, export_path=base_path, export_name_malf=malf_csv, export_name_alf=alf_csv, export_name_mis=days_mis_csv)
malf4
alfs
mis_alfs
min_mis_alfs
min_date_alf

## Estimate
new1, reg = h3.flow_reg(y=sites2, buffer_dis=40000, below_median=True)
reg.to_csv(join(base_path, reg_csv), header=True)

min_flow_malf = new1.malf7d(export_path=base_path, export_name_malf=min_flow_malf_csv)

new2 = new1.data.reset_index().pivot_table(index='time', columns='site', values='data')

##############################################
#### Catchment delineation and naturalisation

catch1 = rec_catch_del(sites_shp=join(base_path, sites_shp), sites_col='site', catch_output=join(base_path, catch_del_shp))

min_flow_nat, min_flow_allo_use = stream_nat(sites2, catch_shp=join(base_path, catch_del_shp), flow_csv=new2, pivot=True, return_data=True, export_path=join(base_path, min_flow_nat_csv))

min_flow_allo_use.to_csv(join(base_path, min_flow_allo_use_csv), header=True, index=False)
min_flow_nat_malf = malf7d(min_flow_nat['nat_flow'], export_path=base_path, export_name_malf=min_flow_nat_malf_csv)

rec_flow_nat, rec_flow_allo_use = stream_nat(sites1, pivot=True, return_data=True, export_path=join(base_path, rec_flow_nat_csv))

rec_flow_allo_use.to_csv(join(base_path, rec_flow_allo_use_csv), header=True, index=False)
rec_nat_malf = malf7d(rec_flow_nat['nat_flow'], export_path=base_path, export_name_malf=rec_nat_malf_csv)

###############################################
#### Plotting

### Flow related

## hydrograph
plt1 = h1.plot_hydrograph(flow_sites=sites1, x_period=x_period, time_format=time_format, start=start, end=end)


#################################################










