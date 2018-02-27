# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas for MALf calcs.
"""

#### Import hydro class
from core.classes.hydro import hydro

#################################################
#### Load data

### Parameters
## Import
mtypes1 = 'flow'
sites1 = [70105, 69607, 69602, 65101, 69505]
qual_codes = [10, 18, 20, 30, 50]

## MALF
intervals = [10, 20]

## Export
base_path = r'S:\Surface Water\temp'
malf_csv = 'malf.csv'
alf_csv = 'alf.csv'
days_mis_csv = 'alf_days_mis.csv'

## Plotting
start = '1986-07-01'
end = '1987-06-30'
x_period = 'month'
time_format = '%d-%m-%Y'

flow_sites = 70105

################################################
#### Import data

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)

################################################
#### Tools

### Flow tools

## MALF and flow stats
fstats = h1.stats(mtypes=mtypes1)
fstats

malf1 = h1.malf7d()
malf1

malf3 = h1.malf7d(intervals=intervals)
malf3

# Get out ALFs and a whole lot of other things
malf4, alfs, mis_alfs, min_mis_alfs, min_date_alf = h1.malf7d(intervals=intervals, return_alfs=True)
malf4
alfs
mis_alfs
min_mis_alfs
min_date_alf

### MALF with csv export
malf4, alfs, mis_alfs, min_mis_alfs, min_date_alf = h1.malf7d(intervals=intervals, return_alfs=True, export=True, export_path=base_path, export_name_malf=malf_csv, export_name_alf=alf_csv, export_name_mis=days_mis_csv)


###############################################
#### Plotting

### Flow related

## hydrograph
plt1 = h1.plot_hydrograph(flow_sites=flow_sites, x_period=x_period, time_format=time_format, start=start, end=end)


#################################################










