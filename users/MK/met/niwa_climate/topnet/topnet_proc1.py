# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:38:19 2017

@author: MichaelEK
"""

from xarray import open_dataset, open_mfdataset, concat
from os.path import join
from core.misc import rd_dir, unarchive_dir
from geopandas import read_file
from core.ecan_io.met import nc_add_gis
from pandas import to_datetime
from core.ts.met import proc_waimak_nc

###########################################
#### Parameters

base_path = r'I:\niwa_data\topnet\waimak2'
out_path = r'E:\ecan\shared\base_data\niwa\climate_projections\topnet\waimak2'

base_path = r'I:\niwa_data\topnet\canterbury'
out_path = r'N:\niwa_netcdf\topnet\canterbury'


##########################################
#### Process data

#unarchive_dir(out_path, 'gz', True)

proc_waimak_nc(base_path, out_path)


























##############################################
#### Testing


#### Parameters

base_path = r'I:\niwa_data\topnet\waimak\RCP4.5'
nc1 = 'streamq_daily_average_1972010200_1985123100_utc_topnet_13070003_strahler3-R1.nc'

nc2 = r'I:\niwa_data\topnet\waimak\RCPpast\BCC-CSM1.1\streamq_daily_average_1986010100_1995123100_utc_topnet_13070003_strahler3-R1.nc'

nc3 = r'I:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1\TotalPrecipCorr_VCSN_BCC-CSM1.1_RCP2.6_2006_2120_south-island_p05_daily_ECan.nc'

nc4 = r'I:\niwa_data\topnet\waimak\RCPpast\CESM1-CAM5\streamq_daily_average_1972010200_1985123100_utc_topnet_13070003_strahler3-F1.nc'

nc5 = r'I:\niwa_data\topnet\waimak\RCP8.5\BCC-CSM1.1\streamq_daily_average_2021010100_2030123100_utc_topnet_13070003_strahler3-S1.nc'

nc6 = r'I:\niwa_data\topnet\waimak\RCP6.0\HadGEM2-ES\streamq_daily_average_2081010100_2090123100_utc_topnet_13070003_strahler3-C1.nc'

nc7 = r'I:\niwa_data\topnet\waimak\RCP8.5\GISS-EL-R\streamq_daily_average_2041010100_2050123100_utc_topnet_13070003_strahler3-P1.nc'

nc8 = r'E:\ecan\shared\base_data\niwa\climate_projections\topnet\waimak\RCP6.0\waimak_topnet_RCP6.0_BCC-CSM1.1.nc'

nc9 = r'E:\ecan\shared\base_data\niwa\climate_projections\topnet\waimak\RCP6.0\waimak_topnet_RCP6.0_CESM1-CAM5.nc'

rec_shp = r'E:\ecan\shared\GIS_base\vector\streams\rec_mfe_cant_no_1st_2nd.shp'
export_shp1 = r'E:\ecan\shared\projects\climate_change\GIS\vector\waimak_rec2.shp'

nc10 = r'I:\niwa_data\topnet\waimak2\RCP2.6\BCC-CSM1.1\streamq_daily_average_2006010200_2010123100_utc_topnet_13070003_strahler3-Q1.nc'


##############################################
#### Read data

rec1 = read_file(rec_shp)

ds1 = open_mfdataset(join(base_path, nc1))
ds1b = open_mfdataset(nc2)
ds2.close()

ds1.close()

rchid = ds1.rchid.to_dataframe().rchid


rec2 = rec1[rec1.NZREACH.isin(rchid)]

rec2.to_file(export_shp1)

ds3 = open_dataset(nc3)
ds3.close()


x_coord = 'longitude'
y_coord = 'latitude'

nc_add_gis(nc3, x_coord, y_coord)



ds1['nrch'] = ds1['rchid']
ds2 = ds1.drop(['basin_lat', 'basin_lon', 'end_lat', 'end_lon', 'rchid', 'time_bnds', 'snwarea_elev_mod', 'snwstor_elev'])
ds3 = ds2.squeeze('nens')
ds3['time'] = to_datetime(to_datetime(ds3['time'].data).date)


ds1b['nrch'] = ds1b['rchid']
ds2b = ds1b.drop(['basin_lat', 'basin_lon', 'end_lat', 'end_lon', 'rchid', 'time_bnds', 'snwarea_elev_mod', 'snwstor_elev'])
ds3b = ds2b.squeeze('nens')
ds3b['time'] = to_datetime(to_datetime(ds3b['time'].data).date)


ds4 = concat([ds3, ds3b], 'time')


flow1 = ds2['river_flow_rate_mod'].copy()

flow2 = flow1.sel(nrch=13036997).to_dataframe()

flow3 = flow1.to_dataframe()

flow3.groupby(level='nrch').sum()['river_flow_rate_mod'].sort_values()

### How to convert date and time to just the date
time1 = flow1.time.data
time2 = to_datetime(to_datetime(time1).date)

flow1['time'] = time2.values


ds7 = open_dataset(nc4)
ds7.close()

ds7 = open_dataset(nc5)
ds7.close()

ds6 = open_dataset(nc6)
ds6.close()

ds7 = open_dataset(nc8)
ds7.close()

ds9 = open_dataset(nc9)
ds9.close()

time1 = to_datetime(ds7.time.data)
time2 = to_datetime(ds9.time.data)


time2[~time2.isin(time1)]


ds10 = open_dataset(nc10)
ds10['nrch'] = ds10['rchid']
ds10.rename({'river_flow_rate_mod': 'flow_rate'}, inplace=True)
ds10a = ds10[['flow_rate', 'drainge', 'ovstn_q', 'soilevp', 'soilh2o', 'aprecip', 'avgtemp', 'zbarh2o', 'canevap', 'potevap']]
tds2a = ds10a.drop(['basin_lat', 'basin_lon', 'end_lat', 'end_lon'])
tds3 = tds2a.squeeze('nens')
time1 = to_datetime(to_datetime(tds3['time'].data).date)
tds3['time'] = time1


ds9.close()





































