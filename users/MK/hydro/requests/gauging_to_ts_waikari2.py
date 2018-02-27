# -*- coding: utf-8 -*-
"""
Example script to estimate a full time series at gauging sites from regresions to recorder sites.
"""

from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.classes.hydro.base import hydro
from core.ts.plot import reg_plot
from os.path import join

############################################
#### Parameters

sites = [65120]

min_years = 10

shp = r'E:\ecan\local\Projects\requests\waikari\zone.shp'

## Export

base_dir = r'E:\ecan\local\Projects\requests\waikari'
reg_csv = '65120_reg.csv'

##################################
#### Read in data

h1 = hydro().get_data('flow_m', sites)
h2 = h1.get_data('flow', shp)

################################
#### Run regression and stats

new1, reg = h2.flow_reg(sites, buffer_dis=40000, min_obs=4, below_median=True)
malf = new1.malf7d()

###############################
#### Export data

reg.to_csv(join(base_dir, reg_csv))
#new_g_malf.to_csv(export_malf)









########################
#### Testing

h0 = hydro().get_data('flow', sites)


































