# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class
from core.classes.hydro import hydro, all_mtypes

#################################################
#### Load data

### Parameters

mtypes1 = 'flow'
mtypes2 = 'flow_m'
mtypes3 = 'swl'
mtypes4 = 'gwl'
mtypes5 = 'gwl_m'
mtypes6 = 'usage'
sites1 = [70105, 69607, 69602, 65101, 69505]
sites2 = [66, 137]
sites3 = ['BT27/5020']
sites4 = ['J38/0774', 'J38/0874', 'J38/0811', 'I39/0033']
qual_codes = [10, 18, 20, 30]
poly = r'S:\Surface Water\backups\MichaelE\Projects\otop\GIS\vector\min_flow\catch1.shp'

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)
h2 = h1.get_data(mtypes=mtypes2, sites=sites2, qual_codes=qual_codes)
h3 = h2.get_data(mtypes=mtypes3, sites=sites1, qual_codes=qual_codes)

gwl1 = hydro().get_data(mtypes=mtypes4, sites=sites3, qual_codes=qual_codes)
gwl2 = hydro().get_data(mtypes=mtypes5, sites=sites3)

use1 = hydro().get_data(mtypes=mtypes6, sites=sites4)

## Find sites based on a polygon shapefile with a 100 m buffer distance (for m_flow)
h4 = hydro().get_data(mtypes=[mtypes1, mtypes2], sites=poly, buffer_dis=100, qual_codes=qual_codes)


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
wy1 = h2.resample(period='water year', n_periods=1, fun='mean')
wy1


#############################################
#### Saving data

export_csv1 = r'S:\Surface Water\temp\test_output1.csv'
export_csv2 = r'S:\Surface Water\temp\test_output2.csv'
export_shp = r'S:\Surface Water\temp\sites_geo.shp'
export_nc = r'S:\Surface Water\temp\test_hydro1.nc'

### Export a time series csv

## Only recorders and pivot
h2.to_csv(export_csv1, mtypes='flow', pivot=True)

## Only gauging data in long format
h2.to_csv(export_csv2, mtypes='flow_m')

### Export a shapfile

## of the gauging sites - combining two functions
h2.sel(mtypes='flow_m').to_shp(export_shp)

## All sites
h2.to_shp(export_shp)

### Make an exact copy of the hydro class object as a netcdf file
h2.to_netcdf(export_nc)


#### Reading the saved data

### csv files
new1 = hydro().rd_csv(export_csv1, time='time', mtypes='flow', dformat='wide')
new1

new2 = hydro().rd_csv(export_csv2, time='time', mtypes='flow_m', sites='site', values='data', dformat='long')
new2

### netcdf files
new3 = hydro().rd_netcdf(export_nc)
new3
new3.geo_loc

###############################################
#### Plotting

### Flow related

## hydrograph
h3 = h2.sel('flow', 70105, start='2000', end='2002')
plt1 = h3.plot_hydrograph(70105, x_period='month', time_format='%Y-%m')

plt2 = h2.plot_hydrograph(x_period='month', x_n_periods=4, time_format='%Y-%m', start='2010', end='2011')

## Regressions
plt2, reg2 = h2.plot_reg('flow', 69602, 'flow_m', 137)
reg2




################################################
#### Fancy stuff ;)

### Use a dictionary of geo attributes to select sites
geo_dict = {'catch_grp': 701, 'cwms': 'Orari-Temuka-Opihi-Pareora'}

h10 = hydro().get_data(mtypes=mtypes1, sites=geo_dict)











