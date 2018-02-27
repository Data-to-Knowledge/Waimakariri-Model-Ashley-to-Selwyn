# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class
from core.classes.hydro import hydro, all_mtypes

#################################################
#### Load data

### Parameters

mtypes1 = 'river_flow_cont_qc'
sites1 = [65104, 64602]
qual_codes = [10, 18, 20, 30, 50]

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)


##################################################
#### Look at the attributes and data contained in the new object

## Look at the general stats of the imported data
print(h1)
print(h2)
h3
h4

## Many other attributes
h2.sites
h2.sites_mtypes
h2.mtypes
h2.mtypes_sites
#h2.site_geo_attr

## The time series data as it is stored in the hydro class object
h2.data

## The geo locations of the sites
h2.geo_loc

#################################################
#### Selections

### Extract time series data as a pandas DataFrame

## Only one site from the flow mtype in long format
df1 = h2.sel_ts(mtypes='flow', sites=69607)
df1

## Both sites from the flow mtype in wide format
df2 = h2.sel_ts(mtypes='flow', pivot=True)
df2

## Select flow sites by a polygon shapefile with a 1 km buffer distance (and pivot)
df3 = h2.sel_ts_by_poly(poly, buffer_dis=1000, mtypes='flow', pivot=True)
df3

## Or more simply by normal brackets
df4 = h2['flow', [69607, 70105]]
df4

### Similar selection, but return a hydro class object instead of a Pandas DataFrame

## Only one site from the flow mtype
h3 = h2.sel(mtypes='flow', sites=69607)
h3

## Two sites from both mtypes
h4 = h2.sel(mtypes=['flow', 'flow_m'], sites=[69607, 66])
h4

## Select flow sites by a polygon shapefile with a 1 km buffer distance
h5 = h2.sel_by_poly(poly, buffer_dis=1000, mtypes='flow')
h5

## Select by geo attributes
#attr_dict = {'catch_grp': 701, 'cwms': 'Orari-Temuka-Opihi-Pareora'}
#h6 = h2.sel_by_geo_attr(attr_dict)
#h6


################################################
#### Tools

### Flow tools

## MALF and flow stats - Applied directly on the original hydro object
fstats = h2.stats(mtypes='flow')
fstats

malf1 = h2.malf7d()
malf1

malf2 = h2.malf7d(sites=[65101, 70105])
malf2

malf3 = h2.malf7d(sites=[65101, 70105], intervals=[10, 20])
malf3

malf4, alfs, mis_alfs, min_mis_alfs, min_date_alf = h2.malf7d(intervals=[10, 20], return_alfs=True)
malf4
alfs
mis_alfs
min_mis_alfs
min_date_alf

## The selection functions can first be applied to the original hydro object
h5 = h2.sel_by_poly(poly, buffer_dis=1000, mtypes='flow')

malf5 = h5.malf7d()
malf5

## Estimate MALFs from gauging sites using nearby recorder sites for regressions
## Two gauging sites with a 40 km buffer distance to determine recorder sites
## Only uses gaugings with the regressions...so far

new1, reg = h2.flow_reg(y=[137, 66], buffer_dis=40000)
malf6 = new1.malf7d()
new1
reg
malf6

## Resampling the time series - returning a hydro class object
wy1 = h1.resample(resample='M', fun='median')
wy1


#############################################
#### Saving data

export_csv1 = r'E:\ecan\local\Projects\requests\jeanine\2017-11-10\64602_65104_mon_median.csv'

### Export a time series csv

## Only recorders and pivot
wy1.to_csv(export_csv1, pivot=True)








