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
sites1 = ['BY20/0018']
#sites2 = ['L36/0633']
qual_codes = [10, 18, 20, 50]
#from_date = '2015-01-01'
#to_date = '2017-06-30'
poly = r'P:\examples\regression_tests\ashburton.shp'

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1)
#h2 = h1.get_data(mtypes=mtypes2, sites=sites2, qual_codes=qual_codes)


## Find sites based on a polygon shapefile with a 10 m buffer distance
h4 = h1.get_data(mtypes=mtypes2, sites=poly, buffer_dis=10, qual_codes=qual_codes)


################################################
#### Tools

### Flow tools

## MALF and flow stats - Applied directly on the original hydro object
gwstats = h4.stats(mtypes='gwl')
gwstats


## Estimate MALFs from gauging sites using nearby recorder sites for regressions
## Two gauging sites with a 40 km buffer distance to determine recorder sites
## Only uses gaugings with the regressions...so far

new1, reg = h4.gwl_reg(y=sites1, buffer_dis=40000, min_yrs=13, min_obs=10)
new1
reg

## Resampling the time series - returning a hydro class object
wy1 = h4.resample(period='water year', n_periods=1, fun='mean')
wy1


#############################################
#### Saving data

base_path = r'P:\examples\regression_tests'

export_csv1 = 'hamishg.csv'
#export_csv2 = 'test_output2.csv'
#export_shp = 'sites_geo.shp'
#export_nc = 'test_hydro1.nc'

### Export a time series csv

## Only recorders and pivot
#h2.to_csv(join(base_path, export_csv1), mtypes='gwl', pivot=True)
#
### Only gauging data in long format
#h2.to_csv(join(base_path, export_csv2), mtypes='gwl_m')
#
#### Export a shapfile
#
### of the gauging sites - combining two functions
#h2.sel(mtypes='gwl_m').to_shp(join(base_path, export_shp))
#
### All sites
#h2.to_shp(join(base_path, export_shp))
#
#### Make an exact copy of the hydro class object as a netcdf file
#h2.to_netcdf(join(base_path, export_nc))

## Only recorders and pivot
new1.to_csv(join(base_path, export_csv1), mtypes='gwl', pivot=True)

#### Reading the saved data

### csv files
#new1 = hydro().rd_csv(join(base_path, export_csv1), time='time', mtypes='flow', dformat='wide')
#new1
#
#new2 = hydro().rd_csv(join(base_path, export_csv2), time='time', mtypes='flow_m', sites='site', values='data', dformat='long')
#new2
#
#### netcdf files
#new3 = hydro().rd_netcdf(join(base_path, export_nc))
#new3
#new3.geo_loc

###############################################
#### Plotting

## Regressions
plt2, reg2 = h4.plot_reg('gwl', 'L37/0023', 'gwl_m', 'BY20/0018')
reg2













