# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class
from core.classes.hydro import hydro, all_mtypes
from core.ts.plot.sw import hydrograph_plot

#################################################
#### Load data

### Parameters

mtypes1 = 'flow'
sites1 = [68002]
qual_codes = [10, 18, 20, 30, 50]
from_date = '2012-01-01'
to_date = '2016-12-31'

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes, from_date=from_date, to_date=to_date)

##################################################
#### plotting

h2 = h1.sel_ts(mtypes='flow', pivot=True)
dayofyear1 = h2.index.dayofyear
year1 = h2.index.year
doy = dayofyear1.to_series().reset_index(drop=True)
year_index = doy.idxmax()
max_year = year1[year_index]
index_new = h2.index[year1 == max_year]

h3 = h2.copy().reset_index()
h3.loc[:, 'year'] = year1
h3.loc[:, 'day'] = dayofyear1

h4 = h3.pivot_table(68002, 'day', 'year')
h4.index = index_new

hydrograph_plot(h4, x_period='month', time_format='%m')


