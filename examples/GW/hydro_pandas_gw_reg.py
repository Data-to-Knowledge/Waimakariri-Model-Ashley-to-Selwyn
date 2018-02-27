# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class and other classes/functions
from core.classes.hydro import hydro, all_mtypes
from os.path import join

#################################################
#### Load data

### Parameters

mtypes1 = 'gwl_m'
mtypes2 = 'gwl'
sites1 = ['K37/3556']
sites2 = ['L36/0633']
qual_codes = [10, 18, 20, 50]
from_date = '2015-01-01'
to_date = '2017-06-30'
poly = r'P:\examples\regression_tests\ashburton.shp'

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=sites2)
h2 = h1.get_data(mtypes=mtypes2, sites=sites2, qual_codes=qual_codes)


## Find sites based on a polygon shapefile with a 10 m buffer distance
h4 = h1.get_data(mtypes=mtypes2, sites=poly, buffer_dis=10, qual_codes=qual_codes)



##################################################
#### Look at the attributes and data contained in the new object

## Look at the general stats of the imported data
print(h1)
print(h2)
h2
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

## Only one site from the mtype in long format
df1 = h2.sel_ts(mtypes='gwl', sites='L36/0633')
df1

## Both sites from the mtype in wide format
df2 = h2.sel_ts(mtypes='gwl', pivot=True)
df2

## Select sites by a polygon shapefile with a 1 km buffer distance (and pivot)
df3 = h2.sel_ts_by_poly(poly, buffer_dis=1000, mtypes='gwl', pivot=True)
df3

## Or more simply by normal brackets
df4 = h2['gwl', ['L36/0633']]
df4

### Similar selection, but return a hydro class object instead of a Pandas DataFrame

## Only one site from the mtype
h3 = h2.sel(mtypes='gwl', sites='L36/0633')
h3

## Two sites from both mtypes
h4 = h2.sel(mtypes=['gwl', 'gwl_m'], sites=['L36/0633', 'K37/3556'])
h4

## Select sites by a polygon shapefile with a 1 km buffer distance
h5 = h2.sel_by_poly(poly, buffer_dis=1000, mtypes='gwl')
h5

## Select by geo attributes
#attr_dict = {'catch_grp': 701, 'cwms': 'Orari-Temuka-Opihi-Pareora'}
#h6 = h2.sel_by_geo_attr(attr_dict)
#h6


################################################
#### Tools

### Flow tools

## MALF and flow stats - Applied directly on the original hydro object
gwstats = h4.stats(mtypes='gwl')
gwstats


## Estimate MALFs from gauging sites using nearby recorder sites for regressions
## Two gauging sites with a 40 km buffer distance to determine recorder sites
## Only uses gaugings with the regressions...so far

new1, reg = h4.gwl_reg(y=sites2, buffer_dis=40000, min_yrs=10)
new1
reg

## Resampling the time series - returning a hydro class object
wy1 = h4.resample(period='water year', n_periods=1, fun='mean')
wy1


#############################################
#### Saving data

base_path = r'P:\examples\regression_tests'

export_csv1 = 'test_output1.csv'
export_csv2 = 'test_output2.csv'
export_shp = 'sites_geo.shp'
export_nc = 'test_hydro1.nc'

### Export a time series csv

## Only recorders and pivot
h2.to_csv(join(base_path, export_csv1), mtypes='gwl', pivot=True)

## Only gauging data in long format
h2.to_csv(join(base_path, export_csv2), mtypes='gwl_m')

### Export a shapfile

## of the gauging sites - combining two functions
h2.sel(mtypes='gwl_m').to_shp(join(base_path, export_shp))

## All sites
h2.to_shp(join(base_path, export_shp))

### Make an exact copy of the hydro class object as a netcdf file
h2.to_netcdf(join(base_path, export_nc))


#### Reading the saved data

### csv files
new1 = hydro().rd_csv(join(base_path, export_csv1), time='time', mtypes='flow', dformat='wide')
new1

new2 = hydro().rd_csv(join(base_path, export_csv2), time='time', mtypes='flow_m', sites='site', values='data', dformat='long')
new2

### netcdf files
new3 = hydro().rd_netcdf(join(base_path, export_nc))
new3
new3.geo_loc

###############################################
#### Plotting

## Regressions
plt2, reg2 = h4.plot_reg('gwl', 'L37/0693', 'gwl_m', 'K37/3556')
reg2













