# -*- coding: utf-8 -*-
"""
Created on Thu May 19 08:55:47 2016

@author: MichaelEK
"""

from numpy import nan, logical_and, where, round, array, in1d, concatenate, isreal, floor, logical_or, sum
from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.ecan_io import rd_sql
from geopandas import read_file, GeoDataFrame

########################################
##### Parameters

#### Fields and names for databases

### Allocation and usage

## Query fields - Be sure to use single quotes for the names!!!

sd_fields = ['Well_No', 'SD1_30', 'SD1_150']

## Equivelant short names for analyses - Use these names!!!

sd_names = ['wap', 'sd_30', 'sd_150']

#### Databases

allo_loc_export_path = r'E:\ecan\shared\GIS_base\vector\allocations\allo_gis.shp'


### Stream depletion

sd_server = "SQL2012PROD05"
sd_database = "Wells"

sd_table = 'Wells.dbo.Well_StreamDepletion_Locations'


#### other parameters

cols1 = ['crc', 'take_type', 'allo_block', 'wap', 'catch_grp', 'status_det', 'use_type', 'geometry']

catch1 = [696]
status1 = ['Issued - Active', 'Issued - s124 Continuance']

export_data = r'E:\ecan\local\Projects\otop\GIS\vector\opihi\opihi_GW_sd.shp'

#export_mon_path = 'C:/ecan/shared/base_data/usage/sd_est_recent_mon_vol.csv'
#export_sd_est_path = 'C:/ecan/shared/base_data/usage/sd_est_all_mon_vol.csv'
#export_reg_path = 'C:/ecan/shared/base_data/usage/sd_est_reg.csv'

########################################
### Read in data

## Allo and usage

allo_gis = read_file(allo_loc_export_path)
allo_gis1 = allo_gis[cols1]

## Stream depletion

sd = rd_sql(sd_server, sd_database, sd_table, sd_fields)
sd.columns = sd_names


#######################################
### Merge tables and prepare data

allo_gis1.loc[:, 'catch_grp'] = to_numeric(allo_gis1.loc[:, 'catch_grp'], errors='coerce')
allo_gis2 = allo_gis1[allo_gis1.catch_grp.isin(catch1)]
allo_gis3 = allo_gis2[allo_gis2.status_det.isin(status1)]
allo_gis3a = allo_gis3[allo_gis3.take_type.isin(['Take Groundwater'])]

allo_gis4 = allo_gis3a.drop_duplicates(['crc', 'take_type', 'allo_block', 'wap', 'catch_grp', 'status_det', 'use_type'])

allo_gis5 = merge(allo_gis4, sd, on='wap', how='left')

## Add in qualitative connection value

qual150 = allo_gis5.sd_150

indexdirect = qual150 >= 90
indexhigh = (qual150 >= 60) & (qual150 < 90)
indexmod = (qual150 >= 40) & (qual150 < 60)
indexlow = (qual150 >= 1) & (qual150 < 40)

allo_gis5['sd_150_qual'] = nan

allo_gis5.loc[indexdirect, 'sd_150_qual'] = 'direct'
allo_gis5.loc[indexhigh, 'sd_150_qual'] = 'high'
allo_gis5.loc[indexmod, 'sd_150_qual'] = 'moderate'
allo_gis5.loc[indexlow, 'sd_150_qual'] = 'low'

qual30 = allo_gis5.sd_30

indexdirect = qual30 >= 90
indexhigh = (qual30 >= 60) & (qual30 < 90)
indexmod = (qual30 >= 40) & (qual30 < 60)
indexlow = (qual30 >= 1) & (qual30 < 40)

allo_gis5['sd_30_qual'] = nan

allo_gis5.loc[indexdirect, 'sd_30_qual'] = 'direct'
allo_gis5.loc[indexhigh, 'sd_30_qual'] = 'high'
allo_gis5.loc[indexmod, 'sd_30_qual'] = 'moderate'
allo_gis5.loc[indexlow, 'sd_30_qual'] = 'low'

## Reorganize data into two sets

allo_sd = allo_gis5.loc[allo_gis5.sd_30_qual.notnull()].reset_index(drop=True)


################################
### Export data

allo_sd.to_file(export_data)





