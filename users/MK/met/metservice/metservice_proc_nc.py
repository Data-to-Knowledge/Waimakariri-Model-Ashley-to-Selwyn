# -*- coding: utf-8 -*-
"""
Created on Fri Jul 07 09:09:11 2017

@author: MichaelEK

Script to process the metservice netcdf files.
"""

from core.ts.met.metservice import proc_metservice_nc, MetS_nc_to_df
from core.spatial import sel_sites_poly, point_interp_ts
from geopandas import read_file
from pandas import to_datetime, concat

##########################################
### Parameters

nc = r'\\fs02\ManagedShares2\Data\Surface Water\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00.nc'
point_shp = r'\\fs02\ManagedShares2\Data\Surface Water\shared\base_data\metservice\testing\Updated_FFC_rf_sites.shp'

point_site_col = 'LAST_HOUR'

time_col = 'time'
y_col = 'y'
x_col = 'x'
data_col = 'precip'
buffer_dis = 16000
digits = 2

rename_dict = {'site': 'MetConnectID', 'model_date': 'PredictionDateTime', 'time': 'ReadingDateTime', 'precip': 'HourlyRainfall'}

export = r'\\fs02\ManagedShares2\Data\Surface Water\shared\base_data\metservice\testing\out_test1.csv'

#########################################
### Process data

## preprocess the nc file to save it as a proper nc file
new_nc = proc_metservice_nc(nc)

## Extract the data from the new nc file
precip, sites, start_date = MetS_nc_to_df(new_nc)

start_date1 = to_datetime(start_date, format='%Y-%m-%d_%H:%M:%S')

## Select the precip data within the buffer area of the points shp
points = read_file(point_shp)
points_poly1 = points.to_crs(sites.crs).buffer(buffer_dis)
sites2 = sel_sites_poly(sites, points_poly1)

precip2 = precip[precip.site.isin(sites2.site)]

## Interpolate the data at the shp points
data = point_interp_ts(precip2, time_col, x_col, y_col, data_col, point_shp, point_site_col, sites.crs, to_crs=None, interp_fun='cubic', agg_ts_fun=None, period=None, digits=2)

########################################
### Output data

## Reformat
data.loc[:, 'model_date'] = start_date1
data1 = data[['site', 'model_date', 'time', 'precip']].copy()
data2 = data1.rename(columns=rename_dict)

##Output
data2.to_csv(export, index=False)
