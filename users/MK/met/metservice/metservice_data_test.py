# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:11:34 2017

@author: michaelek
"""

from xarray import open_dataset, open_mfdataset
from os import path
from core.misc import rd_dir
from numpy import arange, meshgrid, linspace, mgrid, vstack
from netCDF4 import Dataset
import pyproj
from shapely.geometry import Point
from geopandas import GeoDataFrame, GeoSeries
from pycrs.parser import from_epsg_code
from core.misc import select_sites
from pandas import DataFrame, Series, to_datetime


####################################################
#### Load in and look at data

### parameters
base_path = r'E:\ecan\shared\base_data\metservice'

#nc1 = 'nz4kmN_2_2017050918.00.nc'
#nc1 = 'nz4kmN_2_2017051006.00.nc'
nc1 = 'nz8kmN-NCEP_2_2017050918.00.nc'

proj1="+proj=lcc +ellps=GRS80 +lat_1=-60 +lat_2=-30 +lat_0=-60 +lon_0=167.5 +x_0=211921 +y_0=-1221320 +a=6367470 +b=6367470 +no_defs"

out1 = 'test2.shp'

### Read data

lat_coord = 'south_north'
lon_coord = 'west_east'
time_coord = 'Time'
time_var = 'Times'

x1 = open_dataset(path.join(base_path, nc1))

time1 = to_datetime(x1[time_var].data, format='%Y-%m-%d_%H:%M:%S')

nlat = x1.dims[lat_coord]
nlon = x1.dims[lon_coord]
res = int(x1.attrs['DX'])

lat = arange(nlat, dtype='int32') * res
lon = arange(nlon, dtype='int32') * res

x2 = x1.drop(time_var)
x2.coords[time_coord] = ((time_coord), time1)
x2.coords[lat_coord] = ((lat_coord), lat)
x2.coords[lon_coord] = ((lon_coord), lon)

x1.close()


# Problem is that there are no coordinates

d1 = Dataset(nc1, 'r')


###################################################
#### Create x y grid

nlats = 213; nlons = 165
nlats = 519; nlons = 423

# open a the netCDF file for reading.
#ncfile = Dataset(InputFile,'r')
# expected latitudes and longitudes of grid
lats_check = -836.0 + 8.0*arange(nlats,dtype='float32')
lons_check = -211.918 + 8.0*arange(nlons,dtype='float32')
ACPR_check = 900. + arange(nlats*nlons,dtype='float32') # 1d array
ACPR_check.shape = (nlats,nlons) # reshape to 2d array

lat1 = arange(nlats,dtype='int32') * 8000
lon1 = arange(nlons,dtype='int32') * 8000

Y, X = meshgrid(lon1, lat1)
#Y, X = mgrid[0:nlons, 0:nlats]
xy = vstack((X.flatten(), Y.flatten()))

d1 = arange(nlats*nlons, dtype='int32')
d2 = Series(d1, name='site')


geometry = [Point(xy) for xy in zip(xy[1], xy[0])]

gpd = GeoDataFrame(d2, geometry=geometry, crs=proj1)

gpd1 = gpd.to_crs(crs=from_epsg_code(2193).to_proj4())

gpd1.to_file(path.join(base_path, out1))

gpd2 = gpd.to_crs(crs=from_epsg_code(4326).to_proj4())

x_new = gpd2.geometry.apply(lambda s: s.x)
y_new = gpd2.geometry.apply(lambda s: s.y)


#######################################
#### Testing new function

from core.ecan_io import proc_metservice_nc

nc_path = r'E:\ecan\shared\base_data\metservice\nz8kmN-NCEP_2_2017050918.00.nc'
lat_coord = 'south_north'
lon_coord = 'west_east'
time_coord = 'Time'
time_var = 'Times'

proc_metservice_nc(nc_path, lat_coord, lon_coord, time_coord, time_var)


nc1_path = r'E:\ecan\shared\base_data\metservice\nz8kmN-NCEP_2_2017050918.00_corr.nc'

xa1 = open_dataset(nc1_path)

nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051618.00.nc'
lat_coord = 'south_north'
lon_coord = 'west_east'
time_coord = 'Time'
time_var = 'Times'

proc_metservice_nc(nc, lat_coord, lon_coord, time_coord, time_var)


nc1_path = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051618.00_corr.nc'

xa1 = open_dataset(nc1_path)
xa1
xa1.close()

nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00.nc'
nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051612.00.nc'
nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051618.00.nc'

##############################################
#### Extract good nc files and reproject

from xarray import open_dataset
from core.spatial import xy_to_gpd
from pycrs.parser import from_epsg_code
from os import path
from pandas import merge, concat
from numpy import reshape

nc1_path = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051618.00_corr.nc'
base_path = r'E:\ecan\shared\base_data\metservice\testing'
out1 = 'test3.shp'

xa1 = open_dataset(nc1_path)
xa1

df1 = xa1.to_dataframe()
df1a = df1.reset_index(['south_north', 'west_east']).shift(1, freq='H').set_index(['south_north', 'west_east'], append=True)
df0 = concat([df1, df1a], axis=1, join='inner')


df1 = xa1.to_dataframe().reset_index(['south_north', 'west_east'])
df1a = df1.shift(1, freq='H')
df0 = merge(df1.reset_index(), df1a.reset_index(), on=['Time', 'south_north', 'west_east'], how='inner')
df0['precip'] = df0['ACPR_x'] - df0['ACPR_y']
df2 = df0[['Time', 'south_north', 'west_east', 'precip']]
df2.columns = ['time', 'y', 'x', 'precip']

df0.groupby('Time').sum().idxmax()





df2 = df1.groupby('Time')['ACPR'].sum()

time1 = '2017-05-17 07:00:00'
df3 = df2.set_index('time')
df4 = df3.loc[time1]

gpd3 = xy_to_gpd(df4, 'precip', 'x', 'y', proj1)
gpd4 = gpd3.to_crs(crs=from_epsg_code(2193).to_proj4())
gpd4.to_file(path.join(base_path, out1))

xa1.close()





reshape(sites, 2)

Y, X = meshgrid(lon, lat)
#Y, X = mgrid[0:nlons, 0:nlats]
xy = vstack((X.flatten(), Y.flatten())).T


#####################################
#### Save as goetiff




####################################
#### Convert to df and select project area and interpolate

from core.spatial.raster import grid_interp_ts, save_geotiff
from core.ecan_io.met import proc_metservice_nc, ACPR_to_rate, MetS_nc_to_df
from core.spatial import sel_sites_poly, xy_to_gpd
from geopandas import read_file
from rasterio import open as ras_open
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio import transform

nc1_path = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051618.00_corr.nc'
poly = r'S:\Surface Water\backups\MichaelE\Projects\otop\GIS\vector\min_flow\catch1.shp'
poly = r'E:\ecan\shared\base_data\metservice\testing\study_area.shp'

time_col = 'Time'
y_col = 'south_north'
x_col = 'west_east'
data_col = 'precip'
buffer_dis = 10000
grid_res = 500

test4a = r'E:\ecan\shared\base_data\metservice\testing\test5a.shp'
test4b = r'E:\ecan\shared\base_data\metservice\testing\test5b.shp'
test_ras1 = r'E:\ecan\shared\base_data\metservice\testing\test_ras3.tif'
export_path = r'E:\ecan\shared\base_data\metservice\testing\test_ras4.tif'

### Convert ACPR to houly precip
acpr, sites = MetS_nc_to_df(nc1_path)
precip = ACPR_to_rate(acpr)

### Select the sites within the polygon
poly1 = read_file(poly)
sites1 = sites.to_crs(poly1.crs)
sites_sel = sel_sites_poly(sites1, poly, buffer_dis)
sites2 = sites[sites.site.isin(sites_sel.site)]

precip2 = precip[precip.site.isin(sites2.site)]

### Interpolate grid
new_precip = grid_interp_ts(precip2, time_col, x_col, y_col, data_col, grid_res, sites.crs)







### Save one time for comparison
time1 = precip2.groupby('Time')['precip'].mean().idxmax()

out1 = xy_to_gpd(precip2[precip2.Time == time1], data_col, x_col, y_col, sites.crs)
out1 = out1.to_crs(crs=from_epsg_code(2193).to_proj4())
out1.to_file(test4a)

out2 = xy_to_gpd(new_precip[new_precip.time == time1], data_col, 'x', 'y')
out2.to_file(test4b)


x3 = out2.geometry.apply(lambda x: x.x)
y3 = out2.geometry.apply(lambda x: x.y)

trans2 = transform.from_origin(x3.min() - grid_res/2, y3.max() + grid_res/2, grid_res, grid_res)


import numpy as np
x = np.linspace(-4.0, 4.0, 240)
y = np.linspace(-3.0, 3.0, 180)
X, Y = np.meshgrid(x, y)
Z1 = np.exp(-2 * np.log(2) * ((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 1 ** 2)
Z2 = np.exp(-3 * np.log(2) * ((X + 0.5) ** 2 + (Y + 0.5) ** 2) / 2.5 ** 2)
Z = 10.0 * (Z2 - Z1)


z1 = new_precip[new_precip.time == time1].set_index(['y', 'x']).precip
z2 = z1.unstack().values[::-1]

new_dataset = ras_open(test_ras1, 'w', driver='GTiff', height=z2.shape[0], width=z2.shape[1], count=1, dtype=z2.dtype, crs=crs4, transform=trans1)

new_dataset.write(z2, 1)
new_dataset.close()


t1 = ras_open(test_ras1)
crs4 = t1.crs
trans1 = t1.transform

new_dataset = ras_open(test_ras1, 'w', driver='GTiff', height=z2.shape[0], width=z2.shape[1], count=1, dtype=z2.dtype, crs=from_epsg_code(2193).to_proj4(), transform=trans2)
new_dataset.write(z2, 1)
new_dataset.close()


z1 = new_precip.set_index(['y', 'x']).precip
z2 = z1.unstack().values[::-1]


new_dataset = ras_open(test_ras1, 'w', driver='GTiff', height=z2.shape[0], width=z2.shape[1], count=len(time), dtype=z2.dtype, crs=from_epsg_code(2193).to_proj4(), transform=trans2)
for i in range(1, len(time)+1):
    z3 = new_precip[new_precip.time == time[i-1]].set_index(['y', 'x']).precip.unstack().values[::-1]
    new_dataset.write(z3, i)
new_dataset.close()

#### Save geotiff

save_geotiff(new_precip, data_col, 2193, time_col=time_col, export_path=export_path)



##########################################
### More checks...

from core.spatial.raster import grid_interp_ts, save_geotiff
from core.ecan_io.met import proc_metservice_nc, ACPR_to_rate, MetS_nc_to_df
from xarray import open_dataset
from core.spatial import xy_to_gpd
from numpy import arange, meshgrid, vstack
from core.spatial.vector import convert_crs
from pyproj import transform, Proj
from netCDF4 import Dataset
import gdal

nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00.nc'
nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017051612.00.nc'
nc = r'E:\ecan\shared\base_data\metservice\testing\niwa_test1.nc'

nc1 = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00_corr.nc'
nc1 = r'E:\ecan\shared\base_data\metservice\testing\acpcp.2014.nc'
nc2 = r'E:\ecan\shared\base_data\metservice\testing\acpcp_2014a.nc'
nc3 = r'E:\ecan\shared\base_data\metservice\testing\acpcp.mon.ltm.nc'


proj1="+proj=lcc +ellps=GRS80 +lat_1=-60 +lat_2=-30 +lat_0=-60 +lon_0=167.5 +x_0=211921 +y_0=-1221320 +a=6367470 +b=6367470 +no_defs"
proj2 = '+proj=lcc +ellps=GRS80 +lat_1=-60 +lat_2=-30 +lat_0=-60 +lon_0=167.5 +x_0=0 +y_0=0 +a=6367470 +b=6367470 +no_defs'

ds1 = open_dataset(nc1)
ds1
ds1.close()

ds1 = open_dataset(nc3, decode_times=False)
ds1.close()
ds2.close()

gd1 = gdal.Open(nc1)
gd1.GetMetadata()

proc_metservice_nc(nc)

t1 = x2['ACPR'].to_dataframe().reset_index()
t1.groupby('Time').sum()

x1 = arange(165)
y1 = arange(213)
Y, X = meshgrid(lon, lat)
x2 = X.flatten()
y2 = Y.flatten()
xy = vstack((X.flatten(), Y.flatten()))
id1 = arange(165*213)

gpd1 = xy_to_gpd(id1, x2, y2, crs=proj2)

gpd2 = gpd1.to_crs(crs=convert_crs(2193))
convert_crs(4326)

x3, y3 = transform(Proj(proj1), Proj(convert_crs(4326)), x2, y2)

x4 = x3.reshape((213, 165))
y4 = y3.reshape((213, 165))


##############################################
#### Maybe I'll actually finish it this time...

from core.spatial.raster import grid_interp_ts, save_geotiff
from core.ecan_io.met import proc_metservice_nc, MetS_nc_to_df, sel_interp_agg

nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00.nc'
nc1 = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00_corr.nc'
poly = r'E:\ecan\shared\base_data\metservice\testing\study_area.shp'
output_path = r'E:\ecan\shared\base_data\metservice\testing\test_ras5.tif'

time_col = 'time'
y_col = 'y'
x_col = 'x'
data_col = 'precip'
buffer_dis = 10000
grid_res = 500

proc_metservice_nc(nc)

precip, sites = MetS_nc_to_df(nc1)

new_precip = sel_interp_agg(precip, sites, poly, grid_res, data_col, time_col, x_col, y_col, buffer_dis, agg_ts_fun='sum', period='2H', agg_xy=False, output_format='geotiff', output_path=output_path)


##############################################
#### Combine many files

from xarray import open_mfdataset
from os.path import join

base_dir = r'C:\ecan\ftp\metservice\test1'
nc1 = 'nz8kmN-NCEP_2_2017070612_000.00.nc'
nc2 = 'nz8kmN-NCEP_2_2017070612_001.00.nc'

d1 = open_mfdataset(join(base_dir, '*.nc'), concat_dim='times')

d2 = open_mfdataset(join(base_dir, 'nz8kmN-NCEP_2_2017051612.00.nc'))
d3 = open_mfdataset(join(base_dir, nc2))










