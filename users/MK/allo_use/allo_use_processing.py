# -*- coding: utf-8 -*-
"""
Created on Thu May 19 08:55:47 2016

@author: MichaelEK
"""

from numpy import nan, logical_and, where, round, array, in1d, concatenate, isreal, floor, logical_or, sum
from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric, read_hdf
from core.misc import printf, save_df
from core.ecan_io import rd_sql, rd_hydstra_csv
from core.allo_use import allo_proc, allo_ts_proc, allo_errors, allo_use_proc, hist_sd_use, est_use, ros_proc, w_use_proc, allo_gis_proc
from core.ts import w_resample
from core.spatial import pts_sql_join
from geopandas import read_file

########################################
#### Parameters

max_date_allo = '2017-06-30'
max_date_use_est = '2015-06-30'
sql_join_codes = ['swaz_gis', 'catch_gis', 'cwms_gis']
vcn_grid_shp = r'S:\Surface Water\shared\GIS_base\vector\niwa\NIWA_rain_grid_Canterbury.shp'
#vcn_data_path = 'E:/ecan/shared/base_data/precip/VCN_data'

## Export
base_export_path = 'S:/Surface Water/shared/base_data/usage/'

allo_export_path = 'S:/Surface Water/shared/base_data/usage/allo.csv'
allo_loc_export_path = r'S:\Surface Water\shared\GIS_base\vector\allocations\allo_gis.shp'
crc_waps_export_path = 'S:/Surface Water/shared/base_data/usage/crc_waps.csv'
allo_gis_csv = 'S:/Surface Water/shared/base_data/usage/allo_gis.csv'

allo_ts_mon_export = 'S:/Surface Water/shared/base_data/usage/allo_ts_mon_results.h5'
allo_ts_ann_export = 'S:/Surface Water/shared/base_data/usage/allo_ts_ann_results.csv'
allo_use_export = 'S:/Surface Water/shared/base_data/usage/allo_use_ts_mon_results.h5'
use_ros_export_path = 'S:/Surface Water/shared/base_data/usage/allo_use_ros_mon.csv'
ann_use_ros_export_path = 'S:/Surface Water/shared/base_data/usage/allo_use_ros_ann.csv'
usage_mon_est_export_path = 'S:/Surface Water/shared/base_data/usage/allo_est_use_mon.h5'
export_mon_path = 'S:/Surface Water/shared/base_data/usage/sd_est_recent_mon_vol.csv'
export_sd_est_path = 'S:/Surface Water/shared/base_data/usage/sd_est_all_mon_vol.h5'
export_reg_path = 'S:/Surface Water/shared/base_data/usage/sd_est_reg.csv'
export_use_daily = 'S:/Surface Water/shared/base_data/usage/usage_daily_all.h5'
#export_use_mon = 'C:/ecan/shared/base_data/usage/usage_mon.csv'


#######################################
#### Load in data when needed

allo = read_csv(allo_export_path)
allo_ts_mon = read_hdf(allo_ts_mon_export)
#crc_wap = read_csv(crc_waps_export_path)
usage = read_hdf(export_use_daily).reset_index()
allo_use = read_hdf(allo_use_export)
allo_gis = read_csv(allo_gis_csv)
#allo_xy2 = read_file(export_allo_gis)
vcn_grid = read_file(vcn_grid_shp)[['Data VCN_s', 'Network', 'geometry']]
vcn_grid.columns = ['ecan_id', 'niwa_id', 'geometry']
allo_use_ros = read_csv(use_ros_export_path)
usage_est = read_hdf(usage_mon_est_export_path)

########################################
#### Run functions

### allocation processing
allo = allo_proc(export_path=allo_export_path)

### Determine locations based on a single WAP per consent...and stuff...
allo_gis = allo_gis_proc(allo, export_shp=allo_loc_export_path, export_csv=allo_gis_csv)

### Create a monthly/annual time series of allocation volumes
allo_ts_ann = allo_ts_proc(allo_gis, end=max_date_allo, freq='A')
allo_ts_mon = allo_ts_proc(allo_gis, end=max_date_allo, freq='M')

save_df(allo_ts_ann, allo_ts_ann_export, index=False)
save_df(allo_ts_mon, allo_ts_mon_export, index=False)

### Usage
usage = w_use_proc(export=True, export_path=export_use_daily)
allo_use = allo_use_proc(allo_ts_mon, usage.reset_index(), export_path=allo_use_export)

### ROS
allo_use_ros, allo_use_ros_ann = ros_proc(allo_use, export_use_ros_path=use_ros_export_path, export_ann_use_ros_path=ann_use_ros_export_path)

### Estimate usage for missing data from the years 2012 - 2016
usage_est = est_use(allo_use, allo_use_ros, allo_gis, export=True, export_path=usage_mon_est_export_path, date_end=max_date_use_est)

### Spatial joins of all relevant spatial data to the vcsn shp
vcn_grid = read_file(vcn_grid_shp)[['Data VCN_s', 'Network', 'geometry']]
vcn_grid.columns = ['ecan_id', 'niwa_id', 'geometry']
vcn_grid2 = pts_sql_join(vcn_grid, sql_join_codes).dropna(subset=['cwms'])

### Estimate Stream depleting usage for earlier years (1972-2012)
sd_est_all_mon_vol, sd_est_mon_vol, sd_reg = hist_sd_use(usage_est, allo_gis, vcn_grid2, export=True, export_mon_path=export_mon_path, export_sd_est_path=export_sd_est_path, export_reg_path=export_reg_path)


###########################################################
#### Find the errors (that can automatically be found)

#allo_errors(takes, wap, dates, zone, zone_add, takes_names, wap_names, dates_names, zone_names, zone_add_names, irr_names, stock_names, ind_names, pub_names, irr_par, irr_par_names, export_path=base_export_path)

#####################################
### Testing section

sd_all_wap = read_csv(export_sd_est_path)

allo.sort_values('ann_vol', ascending=False)[:30]

allo.sort_values('ann_vol', ascending=False)[30:60]

allo[allo.crc == 'CRC990088.6']
allo[allo.crc == 'CRC950649']

allo[allo.crc == 'CRC021998.1']


takes[takes.crc == 'CRC052033.1']
takes2[takes2.crc == 'CRC052033.1']



allo_ts[allo_ts.crc == 'CRC981180']
CRC981180


CRC905324

allo_ts[allo_ts.crc == 'CRC001135.1']
CRC001135.1



crc_wap_dates[crc_wap_dates.duplicated(subset=['Dates', 'take_type', 'wap'])]
crc_wap_dates[crc_wap_dates.duplicated(subset=['Dates', 'take_type', 'wap'], keep=False)].sort_values(['wap', 'Dates'])
allo_wap1[allo_wap1.duplicated(subset=['Dates', 'take_type', 'wap', 'use_type'], keep=False)].sort_values(['wap', 'Dates'])

crc_wap[crc_wap.duplicated(subset=['take_type', 'wap'], keep=False)].sort_values(['wap'])


allo_wap1[allo_wap1.crc == 'CRC147355']

allo_use1[allo_use1.duplicated(subset=['Dates', 'take_type', 'wap', 'use_type'], keep=False)].sort_values(['wap', 'Dates'])

allo_use2[allo_use2.duplicated(['crc', 'Dates', 'take_type'])]





date1 = '2013-11-30'
crc1 = 'CRC000047.3'
waps1 = ['N34/0142', 'N34/0143', 'N34/0134']


t1 = use1[in1d(use1.wap, waps1)].sort_values('Dates')
tot1 = round(2726.549370 + 9755.096954, 2)

allo_use2[(allo_use2.crc == crc1) & (allo_use2.Dates == date1)]



poly = rd_sql(code='swaz_gis')
poly2 = rd_sql(code='cwms_gis')
join1 = sjoin(allo_xy, poly, how='left', op='within').drop('index_right', axis=1)
join2 = sjoin(join1, poly2, how='left', op='within').drop('index_right', axis=1)


export_sd_hdf = 'C:/ecan/shared/base_data/usage/sd_est_all_mon_vol.h5'
export_sd_nc = 'C:/ecan/shared/base_data/usage/sd_est_all_mon_vol.nc'

from pandas import read_hdf
from xarray import Dataset


t1['swaz'] = t1['swaz'].map(lambda x: x.encode('unicode-escape').decode('utf-8'))


wap1 = 'I39/0033'

w1 = usage[usage.wap == wap1]

w2 = w1.usage
w2.index = w1.date
w2

crc2 = 'CRC093726'
crc_wap[crc_wap.crc == crc2]

crc1 = 'CRC093305'
allo[allo.crc == crc1]

takes[takes.crc == crc1]


### to_date before from_date
allo0 = allo[allo.status_details.isin(['Terminated - Replaced', 'Issued - Active', 'Terminated - Surrendered', 'Terminated - Expired', 'Terminated - Lapsed', 'Issued - s124 Continuance', 'Terminated - Cancelled'])]
allo1 = allo0[['crc', 'from_date', 'to_date']].copy().drop_duplicates()
allo1.loc[:, 'from_date'] = to_datetime(allo1.loc[:, 'from_date'])
allo1.loc[:, 'to_date'] = to_datetime(allo1.loc[:, 'to_date'])

allo1.loc[:, 'diff'] = (allo1.to_date - allo1.from_date).dt.days

a2 = allo1[allo1['diff'] < -30].sort_values('diff', ascending=False)

a2.to_csv('E:/ecan/shared/base_data/usage/bad_fmDate_toDate.csv', index=False)



crc_wap = rd_sql(code='crc_wap_act_acc')


c1 = usage_est.sort_values('crc')

c1[c1.crc == 'CRC000265.2'].wap
allo[allo.crc == 'CRC000345.1']
allo[allo.crc == 'CRC000265.2']

sd[sd.wap == 'J38/0342']




allo_ts_rates = allo_ts_proc(allo_gis, start='2000-01-01', end=max_date_allo, freq='sw_rates')












