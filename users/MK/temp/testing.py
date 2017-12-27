# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 18:16:59 2016

@author: MichaelEK
"""

######################################################
#### 69618 Rockwood correlations to upstream Opihi

from core.classes.hydro import hydro

### Parameters

mtypes1 = 'flow'
mtypes2 = 'flow_m'
y = [69607, 69615, 69616]
x = [69618]
all_sites = [69615, 69616, 69618]

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

h1 = hydro().get_data(mtypes=mtypes1, sites=x)
h2 = h1.get_data(mtypes=mtypes2, sites=y)

### Add both upper sites together

q1 = h2.sel_ts(sites=[69616, 69615], pivot=True)
q2 = (q1[69616] + q1[69615]).dropna()
q2.name = 69617
q3 = DataFrame(q2)
h2 = h2.add_data(q3, mtypes='flow_m', dformat='wide')
h2.get_geo_loc()
y.append(69617)

### Filter dates

h3 = h2.sel(end='1998-12-31')
h4 = h2.sel(start='1999-01-01')

### Regressions

new1, reg1 = h2.flow_reg(y, x, y_mtype='flow_m', x_mtype='flow', min_yrs=1)
new2, reg2 = h3.flow_reg(y, x, y_mtype='flow_m', x_mtype='flow', min_yrs=1)
new3, reg3 = h4.flow_reg(y, x, y_mtype='flow_m', x_mtype='flow', min_yrs=1)


#### Flow comparisons

f1 = hydro().get_data(mtypes=mtypes1, sites=all_sites)

f2 = f1.sel_ts(pivot=True).dropna().interpolate()

### Specific date ranges
start1 = '1999-01-01'
end1 = '1999-06-01'
start2 = '2007-06-01'
end2 = '2007-12-01'


f3 = f2[start1:end1]
f3.plot()

f4 = f2[start2:end2]
f4.plot()


##############################################
#### Hydstra export improvement

from core.ecan_io import rd_sql

### Parameters
server = 'SQL2012PROD03'
db = 'Hydstra'
period_tab = 'PERIOD'
var_tab = 'VARIABLE'
site_tab = 'SITE'
qual_tab = 'QUALITY'

period_cols = ['STATION', 'DATASOURCE', 'VARFROM', 'VARIABLE', 'PERSTART', 'PEREND']
period_names = ['site', 'datasource', 'varfrom', 'varto', 'start', 'end']
var_cols = ['VARNUM', 'VARNAM', 'VARUNIT', 'SHORTNAME']
var_names = ['var_num', 'var_name', 'var_unit', 'var_short_name']
site_cols = ['STATION', 'STNAME', 'SHORTNAME']
site_names = ['site', 'site_name', 'site_short_name']
qual_cols = ['QUALITY', 'TEXT']
qual_names = ['qual_code', 'qual_name']

### Import
period1 = rd_sql(server, db, period_tab, period_cols, where_col='DATASOURCE', where_val=['A'])
period1.columns = period_names

var1 = rd_sql(server, db, var_tab, var_cols)
var1.columns = var_names


##############################################3
#### Hydro testing

from core.classes.hydro import hydro

sites1 = ['BW22/0008', 'BW22/0009', 'L35/0888', 'L35/0716', 'BW22/0070', 'L35/1048', 'L35/0697', 'L35/0565', 'L35/0951', 'L35/0744', 'BW22/0063', 'L35/0721', 'L35/0783', 'L35/1061']

h1 = hydro().get_data(mtypes='usage', sites=sites1)


##############################################
#### Figure out how to make a netcdf GIS worthy

from xarray import open_dataset, Dataset, DataArray
from core.spatial.vector import convert_crs

ex_nc = r'E:\ecan\shared\base_data\metservice\testing\TotalPrecipCorr_VCSN_BCC-CSM1.1_RCPpast_1971_2005_south-island_p05_daily_ECan.nc'

save1 = r'E:\ecan\shared\base_data\metservice\testing\niwa_test1.nc'

ds1 = open_dataset(ex_nc)
ds1
ds1['rain']
ds1['latitude']
ds1['longitude']
ds1['time']

wgs84 = {'grid_mapping_name': "latitude_longitude", 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137, 'inverse_flattening': 298.257223563}

ds2 = ds1.copy()

#time_attrs = ds2.coords['time'].attrs
#time_attrs.update({'units': 'days since 1971-1-1 0:0:0'})
#ds2.coords['time'].attrs = {'units': 'seconds since 1970-1-1 0:0:0', 'long_name': 'time'}
ds3 = ds2.drop(['elevation', 'time_bounds', 'x_index', 'y_index'])
crs4 = {'grid_mapping': 'crs'}
rain_attrs = ds3['rain'].attrs
rain_attrs.update(crs4)
ds3['rain'].attrs = rain_attrs
crs8 = convert_crs(4326, 'netcdf_dict')
ds_crs = DataArray(4326, attrs=crs8, name='crs').to_dataset()

ds5 = ds3.merge(ds_crs)
ds5.to_netcdf(save1)


ds6 = open_dataset(save1)

###############################################
#### Jen usage data request

from core.classes.hydro import hydro

sites4 = ['J38/0774', 'J38/0874', 'J38/0811', 'I39/0033']
mtypes6 = 'usage'
export = r'E:\ecan\local\Projects\requests\jen\2017-05-29\waps_use3.csv'

use1 = hydro().get_data(mtypes=mtypes6, sites=sites4)
use1.to_csv(export, pivot=True)

use2 = use1.sel_ts(sites='I39/0033', pivot=True)
use2.loc['2014-10-30']



##################################################
#### NIWA data extraction and reformating

from xarray import open_dataset, open_mfdataset
from os import path, walk, makedirs
from core.misc import rd_dir
from numpy import in1d
from pandas import read_csv, merge
from core.spatial import xy_to_gpd, sel_sites_poly
from geopandas import read_file
from fnmatch import filter as fil

mtype_name = {'precip': 'TotalPrecipCorr', 'T_max': 'MaxTempCorr', 'T_min': 'MinTempCorr', 'P_atmos': 'MSLP', 'PET': 'PE', 'RH_mean': 'RelHum', 'R_s': 'SurfRad', 'U_z': 'WindSpeed'}
mtype_param = {'precip': 'rain', 'T_max': 'tmax', 'T_min': 'tmin', 'P_atmos': 'mslp', 'PET': 'pe', 'RH_mean': 'rh', 'R_s': 'srad', 'U_z': 'wind'}

vcsn_sites_csv = r'Z:\Data\VirtualClimate\GIS\niwa_vcsn_wgs84.csv'
id_col = 'Network'
x_col = 'deg_x'
y_col = 'deg_y'

param = ['precip', 'PET']
mtype = 'precip'
poly = r'E:\ecan\local\Projects\requests\zeb\2017-06-01\FinalModelDomain.shp'

base_path = r'I:\niwa_data\climate_projections'
out_path = r'Z:\Data\VirtualClimate\NIWA_RCP\climate\waimak'
test_path1 = r'I:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1\TotalPrecipCorr_VCSN_BCC-CSM1.1_RCP2.6_2006_2120_south-island_p05_daily_ECan.nc'
test_path2 = r'I:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1\PE_VCSN_BCC-CSM1.1_RCP2.6_2006_2100_south-island_p05_daily_ECan.nc'

### Test with new functions

from core.ecan_io.met import proc_niwa_rcp, export_rcp_lst

mtypes = ['precip', 'PET']
poly = r'E:\ecan\local\Projects\requests\zeb\2017-06-01\FinalModelDomain.shp'
base_path = r'I:\niwa_data\climate_projections'
export_path = r'Z:\Data\VirtualClimate\NIWA_RCP\climate\waimak'

proc_niwa_rcp(base_path, mtypes, poly, output_fun=export_rcp_lst, export_path=export_path)


##########################
#### Squalarc detection limits

from core.ecan_io import rd_squalarc
from pandas import to_numeric
from numpy import nan

sites = 'E:/ecan/local/Projects/otop/GIS/vector/min_flow/catch1.shp'
sites = ['SQ35872', 'J39/0109']
#sites = ['SQ36017']
from_date = '2010-01-05'
to_date = '2015-06-02'
mtypes = ['Water Temperature']

#t1 = rd_squalarc(sites, mtypes)
t1 = rd_squalarc(sites)
#t1 = rd_squalarc(sites, from_date=from_date, to_date=to_date)
t1 = rd_squalarc(sites, from_date=from_date, to_date=to_date, dtl_method='trend')

########################
#### Improving the MALF function

from core.classes.hydro import hydro
from scipy.interpolate import splrep, splev
from pandas import to_datetime

mtypes1 = 'flow'
sites1 = [70105, 69607, 69602, 65101, 69505]
qual_codes = [10, 18, 20, 30]

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)

df1 = h1.sel_ts(pivot=True)
x = df1.copy()


d1 = df2.loc['1994-07-01':'1995-06-30', 65101]
d2 = df2.loc['1983-07-01':'1984-06-30', 65101]
d3 = df2.loc[:, 65101].values

d4 = df2.interpolate(method='spline', order=2)
d5 = d4.loc['1994-07-01':'1995-06-30', 65101]

dayofyear1 = min_date.apply(lambda x: x.dt.dayofyear)

def day_june(df, dayofyear=182):
    day1 = df.dt.dayofyear
    over1 = day1 >= dayofyear
    under1 = day1 < dayofyear
    day2 = day1.copy()
    day2.loc[over1] = day1.loc[over1] - dayofyear
    day2.loc[under1] = 365 - dayofyear + day1.loc[under1]
    return(day2)

dayofyear1 = min_date.apply(day_june)
dayofyear1.mean()
dayofyear1.std()

series = df2.loc['1994-07-01':'1995-06-30', 65101]

d10 = df2.loc['1975-07-01':'1976-06-30', 65101]
d10 = df2.loc['2008-07-01':'2009-06-30', 65101]
d10 = df2.loc['2001-07-01':'2002-06-30', 70105]
d10 = df2.loc['1986-07-01':'1987-06-30', 70105]

d10.plot()


from core.classes.hydro import hydro

mtypes1 = 'flow'
sites1 = [70105, 69607, 69602, 65101, 69505]
qual_codes = [10, 18, 20, 30]

h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)

malf = h1.malf7d()
malf, alf, alf_mis, alf_min_mis, min_date_alf = h1.malf7d(return_alfs=True)

start = '1986-07-01'
end = '1987-06-30'
flow_sites = 70105

h1.plot_hydrograph(flow_sites=flow_sites, x_period='month', time_format='%d-%m-%Y', start=start, end=end)






###############################################
#### Catchment numbering uniqueness

from geopandas import read_file

shp1 = r'P:\Surface Water Quantity\New surface water catchments\Surface_water_catchments_copy.shp'

g1 = read_file(shp1)

g1.loc[g1['OldHydrolo'].duplicated(), 'OldHydrolo'].sort_values()
g1.loc[g1['OldHydrolo'].isnull(), 'OldHydrolo']


##############################################
#### Graham MALF request

from core.classes.hydro import hydro

mtypes1 = 'flow'
sites3 = [68002, 68001]
qual_codes = [10, 18, 20, 30]

base_path = r'E:\ecan\local\Projects\requests\graham\2017-06-07'
malf_path = '68001_68002_malf.csv'
alf_path = '68001_68002_alf.csv'
alf_mis_path = '68001_68002_alf_missing_days.csv'

h1 = hydro().get_data(mtypes=mtypes1, sites=sites3, qual_codes=qual_codes)
malf, alf, alf_mis, alf_min_mis, min_date_alf = h1.malf7d(return_alfs=True, export_name_malf=malf_path, export_name_alf=alf_path, export_name_mis=alf_mis_path, export_path=base_path, export=True)

d10 = h1.sel_ts(sites=68002, pivot=True, start='1996-07-01', end='1997-06-30')

d10.plot()

x = h1.sel_ts(pivot=True)


###########################################
#### Stony River?

from core.classes.hydro import hydro

mtypes1 = 'flow'
mtypes2 = 'usage'
sites1 = 2280
sites2 = ['I39/0033']

qual_codes = [10, 18, 20, 30]


h1 = hydro().get_data(mtypes=mtypes1, sites=sites1, qual_codes=qual_codes)
h1 = hydro().get_data(mtypes=mtypes2, sites=sites2)


###########################################
#### Marble point

from core.allo_use import restr_days

export1 = r'E:\ecan\local\Projects\requests\suz\2017-06-12\daily_restr_mon.csv'
sites = [64602, 64609, 65104, 65101]

rd1 = restr_days(sites, period='M', export=True, export_path=export1)

##########################################
#### VCSN rearrangement

from xarray import open_dataset, DataArray
from core.spatial.vector import convert_crs

nc = r'Z:\Data\VirtualClimate\vcsn_precip_et_2016-06-06.nc'
x_coord = 'x'
y_coord = 'y'




nc_crs = {'inverse_flattening': 298.25722356299997, 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137, 'transform_name': 'latitude_longitude'}

ds_crs = DataArray(4326, attrs=nc_crs, name='crs').to_dataset()

crs4 = {'grid_mapping': 'crs'}
rain_attrs = ds3['rain'].attrs
rain_attrs.update(crs4)

ds1 = open_dataset(nc)

df1 = ds1.to_dataframe()
df2 = df1.reset_index()

sites_df1 = df2[['x', 'y', 'site']].drop_duplicates().sort_values(['x', 'y']).set_index(['x', 'y'])

data_df1 = df2[['x', 'y', 'time', 'precip', 'ET']].sort_values(['x', 'y', 'time'])
data_df2 = data_df1.set_index(['x', 'y', 'time'])

p1 = data_df2.precip

sites_xa = sites_df1.to_xarray()


##################################################
### Bokeh testing

from bokeh.plotting import figure, output_file, show
from os.path import join

base_dir = r'S:\Surface Water\backups\MichaelE\Projects\plotting\bokeh'
base_dir = r'E:\ecan\local\Projects\plotting\bokeh'
html1 = 'lines.html'
html2 = 'log_lines.html'
html3 = 'color_scatter.html'
html4 = "linked_panning.html"
html5 = "linked_brushing.html"
html6 = "stocks.html"
svg6 = "stocks.svg"
png6 = 'stocks.png'

# prepare some data
x = [1, 2, 3, 4, 5]
y = [6, 7, 2, 4, 5]

# output to static HTML file
output_file(join(base_dir, html1))

# create a new plot with a title and axis labels
p = figure(title="simple line example", x_axis_label='x', y_axis_label='y', logo=None)

# add a line renderer with legend and line thickness
p.line(x, y, legend="Temp.", line_width=2)

# show the results
show(p)


# prepare some data
x = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
y0 = [i**2 for i in x]
y1 = [10**i for i in x]
y2 = [10**(i**2) for i in x]

# output to static HTML file
output_file(join(base_dir, html2))

# create a new plot
p = figure(
   tools="pan,box_zoom,reset,save",
   y_axis_type="log", y_range=[0.001, 10**11], title="log axis example",
   x_axis_label='sections', y_axis_label='particles'
)

# add some renderers
p.line(x, x, legend="y=x")
p.circle(x, x, legend="y=x", fill_color="white", size=8)
p.line(x, y0, legend="y=x^2", line_width=3)
p.line(x, y1, legend="y=10^x", line_color="red")
p.circle(x, y1, legend="y=10^x", fill_color="red", line_color="red", size=6)
p.line(x, y2, legend="y=10^x^2", line_color="orange", line_dash="4 4")

# show the results
show(p)


import numpy as np
from bokeh.plotting import figure, output_file, show

# prepare some data
N = 4000
x = np.random.random(size=N) * 100
y = np.random.random(size=N) * 100
radii = np.random.random(size=N) * 1.5
colors = ["#%02x%02x%02x" % (int(r), int(g), 150) for r, g in zip(50+2*x, 30+2*y)]

# output to static HTML file (with CDN resources)
output_file(join(base_dir, html3), title="color_scatter.py example", mode="cdn")

TOOLS="resize,crosshair,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select"

# create a new plot with the tools above, and explicit ranges
p = figure(tools=TOOLS, x_range=(0,100), y_range=(0,100))

# add a circle renderer with vectorized colors and sizes
p.circle(x,y, radius=radii, fill_color=colors, fill_alpha=0.6, line_color=None)

# show the results
show(p)


########
import numpy as np

from bokeh.layouts import gridplot
from bokeh.plotting import figure, output_file, show

# prepare some data
N = 100
x = np.linspace(0, 4*np.pi, N)
y0 = np.sin(x)
y1 = np.cos(x)
y2 = np.sin(x) + np.cos(x)

# output to static HTML file
output_file(join(base_dir, html4))

# create a new plot
s1 = figure(width=250, plot_height=250, title=None)
s1.circle(x, y0, size=10, color="navy", alpha=0.5)

# NEW: create a new plot and share both ranges
s2 = figure(width=250, height=250, x_range=s1.x_range, y_range=s1.y_range, title=None)
s2.triangle(x, y1, size=10, color="firebrick", alpha=0.5)

# NEW: create a new plot and share only one range
s3 = figure(width=250, height=250, x_range=s1.x_range, title=None)
s3.square(x, y2, size=10, color="olive", alpha=0.5)

# NEW: put the subplots in a gridplot
p = gridplot([[s1, s2, s3]], toolbar_location=None)

# show the results
show(p)


##############
import numpy as np
from bokeh.plotting import *
from bokeh.models import ColumnDataSource

# prepare some date
N = 300
x = np.linspace(0, 4*np.pi, N)
y0 = np.sin(x)
y1 = np.cos(x)

# output to static HTML file
output_file(join(base_dir, html5))

# NEW: create a column data source for the plots to share
source = ColumnDataSource(data=dict(x=x, y0=y0, y1=y1))

TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select,lasso_select"

# create a new plot and add a renderer
left = figure(tools=TOOLS, width=350, height=350, title=None)
left.circle('x', 'y0', source=source)

# create another new plot and add a renderer
right = figure(tools=TOOLS, width=350, height=350, title=None)
right.circle('x', 'y1', source=source)

# put the subplots in a gridplot
p = gridplot([[left, right]])

# show the results
show(p)


###################
import numpy as np

from bokeh.plotting import figure, output_file, show
from bokeh.sampledata.stocks import AAPL
from scipy.signal import fftconvolve
from bokeh.io import export_svgs, export_png

# prepare some data
aapl = np.array(AAPL['adj_close'])
aapl_dates = np.array(AAPL['date'], dtype=np.datetime64)

window_size = 10
window = np.ones(window_size)/float(window_size)
aapl_avg = fftconvolve(aapl, window, 'same')

# output to static HTML file
output_file(join(base_dir, html6), title="stocks.py example")

# create a new plot with a a datetime axis type
p = figure(width=1800, height=800, x_axis_type="datetime")

# add renderers
p.circle(aapl_dates, aapl, size=4, color='darkgrey', alpha=0.2, legend='close')
p.line(aapl_dates, aapl_avg, color='navy', legend='avg')

# NEW: customize by setting attributes
p.title.text = "AAPL One-Month Average"
p.legend.location = "top_left"
p.grid.grid_line_alpha=0
p.xaxis.axis_label = 'Date'
p.yaxis.axis_label = 'Price'
p.ygrid.band_fill_color="olive"
p.ygrid.band_fill_alpha = 0.1

# show the results
show(p)
p.output_backend = 'svg'
export_svgs(p, filename=join(base_dir, svg6))
export_png(p, filename=join(base_dir, png6))


###################################################
#### ftplib testing

from ftplib import FTP
from os import path

### Parameters
dir1 = r'C:\ecan\ftp\metservice\netcdf_combined'

ftp_site = 'ftp1.met.co.nz'
username = 'wrf_rainfall_time_series'
password = 'WrfR@1nF@ll'

### Read local files


### Read ftp files
ftp1 = FTP(ftp_site)
ftp1.login(user=username, passwd=password)
ftp1.cwd('/netcdf_combined')
files1 = ftp1.nlst()

####################################################
### Export rainfall sites for metservice

from core.ecan_io import rd_sql

server = 'SQL2012PROD03'
db = 'MetConnect'
table = 'RainFallPredictionSites'

out_csv = r'E:\ecan\shared\base_data\metservice\metservice_sites.csv'

rsites = rd_sql(server, db, table)
rsites.to_csv(out_csv, index=False)


#################################################
#### Interpolation testing

from core.ts.met.metservice import proc_metservice_nc, MetS_nc_to_df
from core.ts.met.interp import poly_interp_agg
from geopandas import read_file
from pandas import to_datetime, concat
import numpy as np
from scipy.interpolate import griddata
from core.spatial import sel_sites_poly, point_interp_ts


### Parameters

nc = r'E:\ecan\shared\base_data\metservice\testing\nz8kmN-NCEP_2_2017050918.00.nc'
point_shp = r'E:\ecan\shared\base_data\metservice\catch_del\MetConnect_rf_sites.shp'
poly = r'E:\ecan\shared\base_data\metservice\testing\study_area.shp'

time_col = 'time'
y_col = 'y'
x_col = 'x'
data_col = 'precip'
digits = 2
from_crs=None
to_crs=None
interp_fun='cubic'
agg_ts_fun=None
period=None
point_site_col = 'SITENUMBER'
grid_res = 1000
to_crs=2193

test1 = np.array([156914, 5210595])

## preprocess the nc file to save it as a proper nc file
new_nc = proc_metservice_nc(nc)

## extract the data from the new nc file
precip, sites, start_date = MetS_nc_to_df(new_nc)

df = precip.copy()

from_crs = sites.crs

poly1 = read_file(poly).to_crs(from_crs)

sites2 = sel_sites_poly(sites, poly1, 10000)

df = precip.loc[precip.site.isin(sites2.site)]

t = df.groupby('time')['precip'].sum().idxmax()

output1 = grid_interp_ts(df, time_col, x_col, y_col, data_col, grid_res, from_crs, to_crs, interp_fun='cubic', agg_ts_fun=None, period=None, digits=digits)

t1 = output1.groupby('time')['precip'].sum().idxmax()


new_point_df = point_interp_ts(df, time_col, x_col, y_col, data_col, point_shp, point_site_col, from_crs, to_crs, interp_fun, agg_ts_fun=None, period=None, digits=2)



def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2


grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

points = np.random.rand(1000, 2)
values = func(points[:,0], points[:,1])

grid_z2 = griddata(points1, values, (grid_x, grid_y), method='cubic')

import matplotlib.pyplot as plt
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()

x1 = np.array([0.2, 0.4, 0.5, 0.7, 0.7, 0.8])
y1 = np.array([0.3, 0.4, 0.35, 0.35, 0.7, 0.8])

grid_z3 = griddata(points, values, (x1, y1), method='cubic')


xy = np.column_stack((x, y))
z = set1.values

x_new = points.geometry.apply(lambda p: p.x).round(digits).values
y_new = points.geometry.apply(lambda p: p.y).round(digits).values
xy_new = np.column_stack((x_new, y_new))

grid_z4 = griddata(xy, z, xy_new, method='cubic').round(2)


set2 = df2.loc[df2[time_col] == t, ['x', 'y', 'precip']]

gpd5 = xy_to_gpd('precip', 'x', 'y', set2, from_crs)

gpd5.to_file(r'E:\ecan\shared\base_data\metservice\catch_del\test1.shp')

save_geotiff(set2, 'precip', from_crs, x_col='x', y_col='y', time_col=None, export_path=r'E:\ecan\shared\base_data\metservice\catch_del\test2.tif', grid_res=8000)

from time import time

interp1 = Rbf(x, y, z)
z_int = interp1(x_int, y_int)

s1 = time()
z6 = grid_resample(x, y, z, x_int, y_int, digits=2, method='cubic')
e1 = time()
rdf_t = e1 - s1

s1 = time()
z7 = griddata(xy, z, xy_int, method='cubic').round(2)
e1 = time()
griddata_t = e1 - s1

rdf_t/griddata_t



########
### More testing

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl
%matplotlib inline
from scipy.interpolate import (
    LinearNDInterpolator, RectBivariateSpline,
    RegularGridInterpolator, CloughTocher2DInterpolator, SmoothBivariateSpline)
from scipy.ndimage import map_coordinates
#from dolointerpolation import MultilinearInterpolator
# Tweak how images are plotted with imshow
mpl.rcParams['image.interpolation'] = 'none' # no interpolation
mpl.rcParams['image.origin'] = 'lower' # origin at lower left corner
mpl.rcParams['image.cmap'] = 'RdBu_r'

def f_2d(x,y):
    '''a function with 2D input to interpolate on [0,1]'''
    twopi = 2*pi
    return np.exp(-x)*np.cos(x*2*pi)*np.sin(y*2*pi)

def f_3d(x,y,z):
    '''a function with 3D input to interpolate on [0,1]'''
    twopi = 2*pi
    return np.sin(x*2*pi)*np.sin(y*2*pi)*np.sin(z*2*pi)


Ndata = 50
xgrid = np.linspace(0,1, Ndata)
ygrid = np.linspace(0,1, Ndata+1) # use a slighly different size to check differences
zgrid = np.linspace(0,1, Ndata+2)

f_2d_grid = f_2d(xgrid.reshape(-1,1), ygrid)

plt.imshow(f_2d_grid.T)
plt.title(u'image of a 2D function ({}² pts)'.format(Ndata));

f_2d_grid.shape

f_3d_grid = f_3d(xgrid.reshape(-1,1,1), ygrid.reshape(1,-1,1), zgrid)
f_3d_grid.shape

# Define the grid to interpolate on :
Ninterp = 1000
xinterp = np.linspace(0,1, Ninterp)
yinterp = np.linspace(0,1, Ninterp+1) # use a slighly different size to check differences
zinterp = np.linspace(0,1, 5) # small dimension to avoid size explosion


#### Test 1

# Build data for the interpolator
points_x, points_y = np.broadcast_arrays(xgrid.reshape(-1,1), ygrid)
points = np.vstack((points_x.flatten(), points_y.flatten())).T
values = f_2d_grid.flatten()

x_int1, y_int1 = np.meshgrid(xinterp, yinterp)

%timeit f_2d_interp = LinearNDInterpolator(points, values)
f_2d_interp = LinearNDInterpolator(points, values)

%timeit f_2d_interp = CloughTocher2DInterpolator(points, values)
f_2d_interp = CloughTocher2DInterpolator(points, values)

# Evaluate
%timeit f_2d_interp(xinterp.reshape(-1,1), yinterp)

xy_int = np.column_stack((x_int1.flatten(), y_int1.flatten()))




#### Test 2

# Prepare the coordinates to evaluate the array on :
points_x, points_y = np.broadcast_arrays(xinterp.reshape(-1,1), yinterp)
coord = np.vstack((points_x.flatten()*(len(xgrid)-1) , # a weird formula !
                   points_y.flatten()*(len(ygrid)-1)))
coord.shape

f_2d_interp = map_coordinates(f_2d_grid, coord, order=1)


a = np.arange(20.).reshape((5, 4))
a[4][3] = np.nan

map_coordinates(a, [[0.5, 3], [2, 1.99]], order=1, mode='nearest')


#### Test 4
points_x, points_y = np.broadcast_arrays(xgrid.reshape(-1,1), ygrid)
points = np.vstack((points_x.flatten(), points_y.flatten())).T
values = f_2d_grid.flatten()

x_int1, y_int1 = np.meshgrid(xinterp, yinterp)

%timeit spline1 = SmoothBivariateSpline(points_x.flatten(), points_y.flatten(), values, kx=2, ky=2)
spline1 = SmoothBivariateSpline(points_x.flatten(), points_y.flatten(), values, kx=2, ky=2)

spline1([0.5, 1], [0.5, 1])

spline1(x_int1.flatten()[:1000], y_int1.flatten()[:1000])

###########################
### More interp stuff

import numpy as np




#######################################################
#### precip stats testing

from core.classes.hydro import hydro
from pandas import read_csv
from core.ts.met.met_stats import precip_stats

base_dir = r'\\fs02\ManagedShares2\Data\Surface Water\shared\base_data\precip\precip_data.csv'

data = read_csv(base_dir)
data1 = data.set_index(['site', 'time'])

data = precip.sel_ts(mtypes=['precip'])


#################################################
#### Hydstra fixes

from core.ecan_io import rd_sql, rd_hydstra_by_var, write_sql, rd_hydstra_db
from pandas import concat
from os.path import join
from datetime import date, timedelta

sites = [70105, 62103, 62105, 65104, 66429, 66450]
bad_site = [166602, 366425]

for i in sites:
    print(i)
    t1 = rd_hydstra_db(sites=[i])


###############################################
#### Hydrotel mtypes

from core.ecan_io import rd_sql
from pandas import to_datetime, merge, to_numeric, Grouper
from numpy import ndarray, in1d, where
from core.ts.ts import res
from core.misc.misc import time_switch, select_sites

#### Database parameters
server = 'SQL2012PROD05'
database = 'Hydrotel'

objects_tab = 'Hydrotel.dbo.Objects'
objects_col = ['Object', 'Site', 'Name']

objects1 = rd_sql(server, database, objects_tab, objects_col)

t1 = objects1.Name.unique()
t1.sort()
t1.tolist()

c1 = objects1.groupby('Name')['Site'].count()
c1.sort_values()

###############################################
#### Hydrotel all precip sites test

from core.ecan_io import rd_hydrotel

mtype = 'Rainfall'
from_date = '2017-07-21'
to_date = '2017-07-23'

precip1 = rd_hydrotel(mtype=mtype, from_date=from_date, to_date=to_date)

##############################################
#### Relative paths

import patoolib, fnmatch, os


base_dir = r'E:\ecan\git'
test_dir = r'E:\ecan\git\Ecan.Science.Python.Base\core\ts\met\interp.py'

d1 = os.path.relpath(test_dir, base_dir)
path1, file1 = os.path.split(d1)
count1 = len(path1.split(os.path.sep))

# Read in the file
with open(test_dir, 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('core.', ('.' * count1) + 'core.')

# Write the file out again
with open(test_dir, 'w') as file:
  file.write(filedata)



for root, dirs, files in os.walk(folder):
    for filename in fnmatch.filter(files, '*.' + ext):
        print(os.path.join(root, filename))


##############################################
#### Browns rock

from core.classes.hydro import hydro

mtype = 'flow_tel'
site = '66450'

export1 = r'E:\ecan\local\Projects\requests\browns_rock_export\browns_rock_2017-08-03.csv'

b1 = hydro().get_data(mtype, site)
b1.to_csv(export1, pivot=True)

##############################################
#### 2016-2017 usage/allocation ratios

from pandas import read_hdf, read_csv

allo_use_hdf = 'E:/ecan/shared/base_data/usage/allo_use_ts_mon_results.h5'
allo_use_2017_export = r'E:\ecan\local\Projects\requests\Ilja\2017-08-04\allo_use_2016-2017.csv'
allo_csv = 'E:/ecan/shared/base_data/usage/allo.csv'
allo_ts_hdf = 'E:/ecan/shared/base_data/usage/allo_ts_mon_results.h5'
out_data_hdf = r'E:\ecan\local\Projects\hilltop\2017-08_tests\ht_usage_daily.h5'

allo_use1 = read_hdf(allo_use_hdf)
allo_use2 = allo_use1[(allo_use1.date >= '2016-07-01') & (allo_use1.date < '2017-07-01')]
allo_use3 = allo_use2[allo_use2.usage.notnull()]

use1 = allo_use3.groupby(['crc', 'take_type', 'allo_block', 'wap']).sum()
use1['usage-allo_ratio'] = (use1.usage/use1.allo).round(2)
#use1.to_csv(allo_use_2017_export)

t1 = allo_use2.groupby('crc').sum()
t2 = t1[t1.allo > 1000]
t2.usage.notnull().sum()

t3 = t2[t2.usage.notnull()]
t3.allo.sum()/t2.allo.sum()


allo = read_csv(allo_csv)
allo_ts = read_hdf(allo_ts_hdf)
ht_use = read_hdf(out_data_hdf).reset_index()


a1 = use1.reset_index()

crc1 = 'CRC030012'
c1 = a1[a1.crc == crc1]
c1
allo1 = allo[allo.crc == crc1]
allo1.iloc[0]

crc1 = 'CRC972160.2'
c1 = a1[a1.crc == crc1]
c1
allo1 = allo[allo.crc == crc1]
allo1.iloc[0]

crc1 = 'CRC042573.2'
c1 = a1[a1.crc == crc1]
c1
allo1 = allo[allo.crc == crc1]
allo1.iloc[0]

crc1 = 'CRC093599'
c1 = a1[a1.crc == crc1]
c1
allo1 = allo[allo.crc == crc1]
allo1.iloc[0]

a1[a1.wap == c1.wap.values[0]]
allo[allo.wap == c1.wap.values[0]]

c2 = ht_use[ht_use.site == c1.wap.values[0]].set_index('time').drop('site', axis=1)
c2.plot()







##############################################
#### REC network

from geopandas import read_file, GeoDataFrame
import fiona
from core.ecan_io import rd_sql
from time import time
from os.path import join
from shapely.ops import nearest_points
from pandas import concat
from core.spatial.network import find_upstream_rec, extract_rec_catch, agg_rec_catch
from core.spatial.vector import closest_line_to_pts

catch_shp = r'E:\ecan\shared\GIS_base\vector\catchments\river-environment-classification-watershed-canterbury-2010.shp'
streams_shp = r'E:\ecan\shared\GIS_base\vector\streams\rec-canterbury-2010.shp'

server = 'SQL2012PROD05'
db = 'GIS'
table = 'MFE_NZTM_RECWATERSHEDCANTERBURY'
cols = ['NZREACH']

base_dir = r'P:\cant_catch_delin\set2'
sites_shp = 'sites.shp'
sites_col = 'site'
catch_out = 'catch1.shp'

streams1 = read_file(streams_shp)
catch1 = read_file(catch_shp)

pts1 = read_file(join(base_dir, sites_shp))

pts_seg = closest_line_to_pts(pts1, streams1, line_site_col='NZREACH', dis=1000)
nzreach = pts_seg.copy().NZREACH.unique()

reaches = find_upstream_rec(nzreach)

rec_catch = extract_rec_catch(reaches)

rec_shed = agg_rec_catch(rec_catch)
rec_shed.columns = ['NZREACH', 'geometry']
rec_shed1 = rec_shed.merge(pts_seg[['site', 'NZREACH']], on='NZREACH')

rec_shed1.to_file(join(base_dir, catch_out))

t1 = rec_shed.loc[[13059458]]


catch2 = catch1[catch1.NZREACH.isin(sites)].dissolve('NZREACH')[['geometry']]


##################################################
#### Lake
streams_shp = r'E:\ecan\shared\GIS_base\vector\streams\river-environment-classification-canterbury-2010.shp'
sites_pts_shp = r'P:\Surface Water Quantity\Projects\Freshwater Report\lake1.shp'
catch_shp = r'P:\Surface Water Quantity\Projects\Freshwater Report\lake_catch.shp'

streams = read_file(streams_shp)
pts = read_file(sites_pts_shp)

pts_seg = closest_line_to_pts(pts, streams, line_site_col='NZREACH', dis=50)
nzreach = pts_seg.copy().NZREACH.unique()
#pts_seg.to_file(join(base_dir, rec_sites_rec_shp))

reaches = find_upstream_rec(nzreach)

rec_catch = extract_rec_catch(reaches)
#rec_catch.to_file(join(base_dir, catch_del_temp_shp))

rec_shed = agg_rec_catch(rec_catch)
rec_shed.columns = ['NZREACH', 'geometry', 'area']
rec_shed1 = rec_shed.merge(pts_seg[['id', 'NZREACH']], on='NZREACH')

rec_shed1.to_file(catch_shp)
rec_shed1.area.sum()


############################################
### Manual GW readings

from core.classes.hydro import hydro

wap = 'BT27/5020'

gw1 = hydro().get_data('gwl_m', sites=wap)



##########################################
### New Hilltop connection tests

import Hilltop
from pandas import DataFrame, concat, to_datetime
from core.ecan_io.hilltop import ht_sites, ht_get_data

hts1 = r'H:\Data\Archive\Boraman2016-17.hts'
hts2 = r'H:\Data\Telemetry\Boraman.hts'
hts = r'H:\Data\Squalarc.hts'

csv1 = r'E:\ecan\local\Projects\requests\squalarc\squalarc_site_info.csv'

sites = ['M35/0132', 'SQ32943']

dfile1 = Hilltop.Connect(hts1)
sites = Hilltop.SiteList(dfile1)

site_info = DataFrame()

for i in sites:
    try:
        info1 = Hilltop.MeasurementList(dfile1, i)
    except SystemError:
        print('Site ' + str(i) + " didn't work")
    info1.loc[:, 'site'] = i
    site_info = concat([site_info, info1])
site_info.reset_index(drop=True, inplace=True)

site_info.loc[:, 'Start Time'] = to_datetime(site_info.loc[:, 'Start Time'], format='%d-%b-%Y %H:%M:%S')
site_info.loc[:, 'End Time'] = to_datetime(site_info.loc[:, 'End Time'], format='%d-%b-%Y %H:%M:%S')

index1 = 1

d1 = Hilltop.GetData(dfile1, site_info.loc[index1, 'site'], site_info.loc[index1, 'Measurement'], '', '', method='Average', interval='1 day', alignment = '00:00')

d1 = Hilltop.GetData(dfile1, site_info.loc[index1, 'site'], site_info.loc[index1, 'Measurement'], '', '')


s1 = ht_sites(hts)
s1.to_csv(csv1, index=False, header=True)

s2 = ht_sites(hts2)

data, miss1 = ht_get_data(hts, agg_method='', interval='', alignment='', output_missing_sites=True)

Hilltop.Disconnect(dfile1)

###############################################
#### Lake levels

from core.ecan_io import rd_hydrotel
from core.classes.hydro import hydro
from core.ecan_io import rd_hydstra_db

site = [71178]

#l1 = rd_hydrotel(site, mtype='swl_tel')


hy1 = hydro()._rd_hydstra(site, datasource='COMBINED', varfrom=140)

hy2 = rd_hydstra_db(sites=site, datasource='FLOW', varfrom=140)
hy2 = rd_hydstra_db(sites=[70105], datasource='A')

with hyd as h:
    var1 = h.get_variable_list(site, 'COMBINED')



############################################
#### 7-day Volume estimate

from core.classes.hydro import hydro
from pandas import concat, TimeGrouper


site = [70105]
mtype = 'flow'

h1 = hydro().get_data(mtype, site)
h2 = h1.sel_ts('flow', pivot=True)

h3 = h2.rolling(7).sum()
h4 = h3.rolling(7, min_periods=1).max()

h3.columns = ['sum']
h4.columns = ['max']

h5 = concat([h2, h3, h4], axis=1)
h5['limit'] = 100

h5['diff'] = (h5['sum'] - h5['limit'])
h5.loc[h5['diff'] < 0, 'diff'] = 0

mon_grp = h5['diff'].rolling(14)


def take_max(arr):
    max1 = arr.max()
    end_val = arr[-1]
    if end_val == max1:
        val = max1
    else:
        val = 0
    return(val)


mon_grp.apply(take_max)


xy1 = precip_et[precip_et.time == precip_et.time[0]]

xy2 = xy_to_gpd(xy1.index, 'x', 'y', xy1, 4326)
xy2.columns = ['site', 'geometry']
xy2.to_file(r'D:\lsrm_results\test\test5.shp')

from pandas import read_hdf

t1 = r'D:\lsrm_results\RCP4.5_BCC-CSM1.1_80perc.h5'
t2 = read_hdf(t1)


##############################################
#### Squalarc dup parameters

param1 = rd_sql('SQL2012PROD05', 'Squalarc', '"SQL_SAMPLE_METHODS+"', stmt='select distinct PA_NAME from Squalarc.dbo."SQL_SAMPLE_METHODS+"')
param1.columns = ['param']
param2 = param1.param

param2[param2.str.contains('nitrogen', case=False, na=False)].sort_values()

param2[param2.str.contains('nitrate', case=False, na=False)].sort_values()
param2[param2.str.contains('nitrite', case=False, na=False)].sort_values()


#############################################
#### Use SQLalchemy -- Not yet...

from sqlalchemy import create_engine

stmt = "SELECT * FROM F_WQ_SampleResults"
server = "SQL2012PROD03"
database = "DataWarehouse"

eng = create_engine('mssql+pyodbc://' + server + '/' + database)

##############################################
#### Hydstra list modified

add_where = "TIDEDA_FLAG='N'"
data_col= 'DEPTH_TO_WATER'
database = 'Wells'
qual_col = None
server = 'SQL2012PROD05'
site_col = 'WELL_NO'
table = 'DTW_READINGS'
time_col = 'DATE_READ'
sites = ['BY20/0018']
mtype = 'gwl_m'

#############################################
## New mssql time series agg

from core.ecan_io.mssql import rd_sql_ts

resample = 'hour'
period = 4


stmt1 = "SELECT Point AS site, DT AS time, SampleValue AS value FROM " + data_tab + " where " + " and ".join(where_lst)

data1 = rd_sql(server, database, data_tab, stmt=stmt1).drop('site', axis=1)
data1.set_index(['time'], inplace=True)

day1 = data1.resample('D').mean().round(3)


stmt1 = "SELECT " + "Point, DATEADD(" + resample + ", DATEDIFF(" + resample + ", 0, DT)/ " + str(period) + " * " + str(period) + ", 0) AS time, round(" + fun + "(SampleValue), 3) AS value" + " FROM " + data_tab + " where " + " and ".join(where_lst) + " GROUP BY Point, DATEADD(" + resample + ", DATEDIFF(" + resample + ", 0, DT)/ " + str(period) + " * " + str(period) + ", 0) ORDER BY Point, time"

data2 = rd_sql(server, database, data_tab, stmt=stmt1).drop('site', axis=1)
data2.set_index(['time'], inplace=True)

data2

server = 'SQL2012PROD05'
database = 'Hydrotel'
groupby_cols = 'Point'
date_col = 'DT'
values_cols = 'SampleValue'
resample_code = 'D'
period = 1
val_round=3
fun = 'sum'
table = 'Samples'
from_date = '2017-01-01'
to_date = '2017-09-01'
where_val=None
where_op='AND'
where_col = {'Point': [3333]}

df1 = rd_sql_ts(server, database, table, groupby_cols, date_col, values_cols, resample_code, period, fun, val_round, where_col, where_val, where_op, from_date, to_date)


#################################################
##
from core.ecan_io import rd_hydrotel

py_dir = r'E:\ecan\git\Ecan.Science.Python.Base\tests\classes\hydro'

sites = ['L36/0633']
sites = ['K38/1081']
mtype = 'aq_wl_cont_raw'
from_date=None
to_date=None
resample_code='D'
period=1
fun='mean'
val_round=3

t1 = rd_hydrotel(sites, 'aq_wl_cont_raw')

map1 = {'flow': 'river_flow_cont_qc', 'swl': 'river_wl_cont_qc', 'flow_m': 'river_flow_disc_qc'}

m2 = [map1[i] for i in m1]

qual_codes = [10, 18, 20, 30, 50]
from_date = '2015-07-01'
to_date = '2016-06-30'

h2 = hydro().get_data(mtypes='river_flow_cont_qc', sites=[65101, 69505, 69602, 69607, 70105], qual_codes=qual_codes, from_date=from_date, to_date=to_date)
h2 = h2.get_data(mtypes='river_flow_disc_qc', sites=[66, 137], from_date=from_date, to_date=to_date)
h2 = h2.get_data(mtypes='river_wl_cont_qc', sites=[65101, 69505, 69602, 69607, 70105], qual_codes=qual_codes, from_date=from_date, to_date=to_date)
h2.to_netcdf(path.join(py_dir, netcdf1))

site = 'O31/0372'
time1 = '2005-02-17'



#####################################################
### All lowflows sites

from core.allo_use import crc_band_flow

csv = r'E:\ecan\local\Projects\requests\jeanine\2017-11-09\LF_sites.csv'

site_lf = crc_band_flow()
site_lf.to_csv(csv, index=False)


####################################################
### hydr import issues

from core.classes.hydro import hydro


sites = [1688218, 70105]

server = 'SQL2012PROD03'
database = 'DataWarehouse'
table = 'F_HY_SWL_data'
mtype = 'river_wl_cont_qc'
date_col = 'time'
site_col = 'site'
data_col = 'data'
qual_col='qual_code'
from_date='2016-01-01'
to_date='2017-01-01'
qual_codes=[10, 18, 20, 30, 50]
add_where=None
min_count=20
resample_code='W'
period=1
fun='mean'
val_round=3
where_col = {'qual_code': [10, 18, 20, 30, 50], 'site': [1688218, 70105]}

t1 = rd_hydrotel(sites, from_date=from_date, to_date=to_date, min_count=368)

rd_sql_ts(server, database, table, site_col, date_col, data_col, resample_code, period, fun, val_round, where_col, from_date, to_date, min_count)

i = 'atmos_precip_cont_qc'


sites_sql_fun=geo_loc_dict[i]
mtype_dict=mtypes_sql_dict[i]
mtype='atmos_precip_cont_qc'
sites=[313710]
from_date = '2017-02-01'
to_date = '2017-04-01'
qual_codes = [10, 18, 20, 30, 50, 11, 21]
min_count=20
buffer_dis=0
resample_code=None
period=1
fun='mean'

h1 = hydro()

sites=sites3
mtype=mtype
from_date=from_date
to_date=to_date
qual_codes=qual_codes
min_count=min_count
resample_code=resample_code
period=period
fun=fun

**mtype_dict

add_where =  None
data_col =  'data'
database =  'DataWarehouse'
date_col =  'time'
qual_col =  'qual_code'
server =  'SQL2012PROD03'
site_col =  'site'
table =  'F_HY_Precip_data'

database= 'DataWarehouse'
date_col ='time'
from_date ='2017-02-01'
fun ='mean'
groupby_cols= 'site'
min_count =20
period =1
resample_code= None
server= 'SQL2012PROD03'
table ='F_HY_Precip_data'
to_date= '2017-04-01'
val_round= 3
values_cols= 'data'
where_col ={'qual_code': [10, 18, 20, 30, 50, 11, 21], 'site': ['313710']}


h1 = hydro().get_data(mtypes='river_wl_cont_qc', sites=sites, qual_codes=[10, 18, 20, 30, 50, 11, 21], from_date = '2016-02-01', to_date = '2017-04-01', min_count=52, resample_code='W')

##################################################
### Hydrotel

t1 = rd_hydrotel([66204], mtype='river_flow_cont_raw', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'))

##################################################
### Add MetConnect sites table

from pandas import read_csv, to_datetime
from core.ecan_io import write_sql

server = 'SQL2012DEV01'
database = 'Hydro'
database = 'MetConnect'
table = 'RainFallPredictionSitesGrid'

dtype_dict = {'MetConnectID': 'INT', 'SiteString': 'VARCHAR(34)', 'Office': 'VARCHAR(2)', 'HydroTelPointNo': 'INT', 'TidedaID': 'INT', 'StartDate': 'DATE'}

sites_csv = r'E:\ecan\shared\projects\metservice_processing\RainFallPredictionSites.csv'

sites = read_csv(sites_csv)
sites['StartDate'] = to_datetime(sites['StartDate'], dayfirst=True, format='%d/%m/%Y')

stmt_dict = write_sql(sites, server, database, table, dtype_dict, primary_keys=['MetConnectID'])


####################################################
### Low flows

siteid = 12

restr_val[restr_val.SiteID == siteid]
sites[sites.SiteID == siteid]

####################################################
### More hydstra
from pint import UnitRegistry
from core.ecan_io.hydllp import rd_hydstra_by_mtype
from core.ecan_io import rd_sql
from pandas import to_numeric

site = 164606

from_date = '2017-10-01'
to_date = '2017-12-06'

period2.varto.sort_values().unique().tolist()

(period2.varto == 500).sum()
period2[(period2.varto == 101)]

ureg = UnitRegistry()

m_s = ureg.meter / ureg.second

mtype = 'lake_wl_rec_qc'
mtype = 'river_flow_rec_qc'

from_date = '2017-01-01'
to_date = '2017-01-10'
interval = 'hour'
data_type = 'point'


data1 = rd_hydstra_by_mtype(mtype, from_date=from_date, interval=interval)

server= 'SQL2012PROD03'
database = 'Hydstra'
table = 'RATEPER'
fields = ['STATION', 'VARFROM', 'VARTO', 'SDATE', 'STIME', 'REFSTN', 'REFTAB', 'PHASE']

sites = [66612]
from_date = '2017-09-01'
to_date = '2017-12-01'



rate1 = rd_sql(server, database, table, fields)
rate1['STATION'] = rate1['STATION'].str.strip()
rate1['REFSTN'] = rate1['REFSTN'].str.strip()
rate1['STATION'] = to_numeric(rate1['STATION'], 'coerce')
rate1['REFSTN'] = to_numeric(rate1['REFSTN'], 'coerce')

rate2 = rate1[rate1.STATION == site]
rate2.sort_values('SDATE')

#########################################
### pint tests

from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity












