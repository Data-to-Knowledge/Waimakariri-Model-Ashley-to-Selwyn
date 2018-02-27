# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 16:27:09 2017

@author: MichaelEK
"""

import sys
sys.path.append(r'C:\git\Ecan.Science.Python.Base')

from geopandas import read_file, sjoin, GeoDataFrame
from pandas import DateOffset, to_datetime, concat, merge, cut, DataFrame, MultiIndex, Series, read_csv
from core.spatial.vector import multipoly_to_poly, spatial_overlays
from core.ts import grp_ts_agg
from datetime import date
from scipy.stats import percentileofscore
from numpy import nan
from core.classes.hydro import hydro
from core.ts import tsreg
from core.ecan_io import rd_hydrotel
from configparser import ConfigParser
from os import path, getcwd
from warnings import filterwarnings

from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, CategoricalColorMapper, CustomJS, renderers, annotations, Panel, Tabs
from bokeh.palettes import brewer
from bokeh.models.widgets import Select
from bokeh.layouts import column

from Tkinter import Tk
from tkFileDialog import askdirectory

date1 = date.today()
date_str = str(date1)
filterwarnings('ignore')
Tk().withdraw()

###################################################
#### Parameters

base_dir = r'P:\Surface Water Quantity\Projects\Freshwater Report\sw_precip'

sw_poly_shp = 'sw_boundary_v01.shp'
precip_poly_shp = 'precip_boundary_v01.shp'
#gw_poly_shp = 'precip_boundary_v01.shp'
rec_catch_shp = 'catch_delin_recorders.shp'
#gw_sites_shp = 'gw_sites1.shp'
view_bound_shp = 'boundary_v02.shp'
precip_site_shp = 'precip_sites.shp'
pot_sw_site_list_csv = 'potential_sw_site_list.csv'

qual_codes = [10, 18, 20, 30, 50, 11, 21, 40]

n_previous_months = 6

month_names = ['Jan', 'Feb', 'March', 'April', 'May', 'June', 'July', 'August', 'Sept', 'Oct', 'Nov', 'Dec']

lon_zone_names = {'L': 'Lowlands', 'F': 'Foothills', 'M': 'Mountains', 'BP': 'Banks Peninsula'}

### Output
output_dir = askdirectory(initialdir=base_dir, title='Select the output directory', mustexist=True)
ts_out_csv = 'ts_out_perc_' + date_str + '.csv'

## plots
#precip_sw_gw1_html = 'precip_sw_gw_' + date_str + '.html'
precip_sw1_html = 'precip_sw_' + date_str + '.html'
sw1_html = 'sw_catch_' + date_str + '.html'

##################################################
#### Read in data

print('Reading in the data')

### Overall veiw
view_zones = read_file(path.join(base_dir, view_bound_shp))
view_zones = view_zones.replace({'lon_zone': lon_zone_names})
view_zones['zone'] = view_zones['lat_zone'] + ' - ' + view_zones['lon_zone']
view_zones = view_zones.drop(['lon_zone', 'lat_zone'], axis=1)

### SW
sw_poly = read_file(path.join(base_dir, sw_poly_shp))[['lat_zone', 'lon_zone', 'geometry']]
rec_catch = read_file(path.join(base_dir, rec_catch_shp))
site_list = read_csv(path.join(base_dir, pot_sw_site_list_csv))

sw_list = site_list.replace({'lon_zone': lon_zone_names})
sw_list['zone'] = sw_list['lat_zone'] + ' - ' + sw_list['lon_zone']
sw_list = sw_list.drop(['lon_zone', 'lat_zone'], axis=1)
sw_list.loc[19:20, 'Notes'] = nan

sw_zones = sw_poly.replace({'lon_zone': lon_zone_names})
sw_zones['zone'] = sw_zones['lat_zone'] + ' - ' + sw_zones['lon_zone']
sw_zones = sw_zones.drop(['lon_zone', 'lat_zone'], axis=1)
sw_zones['mtype'] = 'flow'

### precip
precip_sites = read_file(path.join(base_dir, precip_site_shp))
precip_zones = read_file(path.join(base_dir, precip_poly_shp))

precip_zones = precip_zones.replace({'lon_zone': lon_zone_names})
precip_zones['zone'] = precip_zones['lat_zone'] + ' - ' + precip_zones['lon_zone']
precip_zones = precip_zones.drop(['lon_zone', 'lat_zone'], axis=1)
precip_zones['mtype'] = 'precip'

### gw
#gw_sites = read_file(path.join(base_dir, gw_sites_shp))
#gw_zones = read_file(path.join(base_dir, gw_poly_shp))
#
#gw_zones = gw_zones.replace({'lon_zone': lon_zone_names})
#gw_zones['zone'] = gw_zones['lat_zone'] + ' - ' + gw_zones['lon_zone']
#gw_zones = gw_zones.drop(['lon_zone', 'lat_zone'], axis=1)
#gw_zones['mtype'] = 'gw'
#
#gw_sites = gw_sites.replace({'lon_zone': lon_zone_names})
#gw_sites['zone'] = gw_sites['lat_zone'] + ' - ' + gw_sites['lon_zone']
#gw_sites = gw_sites.drop(['lon_zone', 'lat_zone'], axis=1)


### Combine

#zones = concat([sw_zones, precip_zones, gw_zones]).reset_index(drop=True)
zones = concat([sw_zones, precip_zones]).reset_index(drop=True)

#################################################
#### Select sites

### SW
sites1 = sw_list[sw_list.Notes.isnull()].drop('Notes', axis=1)

flow1 = hydro().get_data('river_flow_cont_qc', sites1.site, qual_codes)
stats_flow = flow1.stats('river_flow_cont_qc')

### precip
precip1 = hydro().get_data(mtypes='atmos_precip_cont_qc', sites=precip_sites.site, qual_codes=qual_codes)

### GW
#gw1 = hydro().get_data(mtypes='gwl', sites=gw_sites.site, qual_codes=qual_codes)

#################################################
#### Run monthly summary stats

print('Processing past data')

### SW
flow2 = flow1.sel_ts(mtypes='river_flow_cont_qc')
flow2.index = flow2.index.droplevel('mtype')
flow3 = flow2.reset_index()
flow3 = flow3.replace({'site': {66403: 66442}})

mon_flow1 = grp_ts_agg(flow3, 'site', 'time', 'M').median().reset_index()
mon_flow1['mon'] = mon_flow1.time.dt.month
mon_flow1['mtype'] = 'flow'

### precip
precip2 = precip1.sel_ts(mtypes='atmos_precip_cont_qc')
precip2.index = precip2.index.droplevel('mtype')
precip3 = precip2.reset_index()

mon_precip1 = grp_ts_agg(precip3, 'site', 'time', 'M').sum().reset_index()
mon_precip1['mon'] = mon_precip1.time.dt.month
mon_precip1['mtype'] = 'precip'

### gw
#gw2 = gw1.sel_ts(mtypes='gwl')
#gw2.index = gw2.index.droplevel('mtype')
#gw3 = gw2.reset_index()
#
#mon_gw1 = grp_ts_agg(gw3, 'site', 'time', 'M').median().reset_index()
#mon_gw1['mon'] = mon_gw1.time.dt.month
#mon_gw1['mtype'] = 'gw'

### Combine all mtypes

#mon_summ = concat([mon_flow1, mon_precip1, mon_gw1]).reset_index(drop=True)
mon_summ = concat([mon_flow1, mon_precip1]).reset_index(drop=True)

###############################################
#### Pull out recent monthly data from hydrotel

now1 = to_datetime(date1)
start_date = now1 - DateOffset(months=n_previous_months) - DateOffset(days=now1.day - 1)
end_date = now1 - DateOffset(days=now1.day - 1)

### SW
print('Getting HydroTel Flow Data:')
sites2 = sites1.copy()
sites2.loc[sites2.site.isin([64610, 65104, 68526]), 'site'] = [164610, 165104, 168526]

hy_sites = sites2.site.copy()

hy1 = rd_hydrotel(hy_sites, mtype='river_flow_cont_raw', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample_code='D', fun='mean')
hy2 = hy1.reset_index()
if len(hy2.site.unique()) != len(hy_sites):
    print(str(len(hy_sites) - len(hy2.site.unique())) + " sites are not in Hydrotel")
last1 = hy2.groupby('site').last()
last_index = (last1.time >= end_date - DateOffset(days=1))
if not all(last_index):
    print(str(sum(~last_index)) + " sites have less than a full months record")
    hy2 = hy2[hy2.site.isin(last_index.index[last_index])]
hy2 = hy2[hy2.time != end_date]
hy2 = hy2.replace({'site': {66403: 66442}})
hy3 = grp_ts_agg(hy2, 'site', 'time', 'M').median().reset_index()
#hy3.columns = ['site', 'mon_median_flow']

hy3.loc[:, 'site'] = hy3.site.replace([164610, 165104, 168526], [64610, 65104, 68526])

hy_flow = hy3.copy()
hy_flow['mtype'] = 'flow'

### precip
print('Getting HydroTel Precip Data:')

hy_sites = mon_precip1.site.unique().copy()
hy1 = rd_hydrotel(hy_sites, mtype='atmos_precip_cont_raw', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample_code='D', fun='sum')
hy2 = hy1.reset_index()
if len(hy2.site.unique()) != len(hy_sites):
    print(str(len(hy_sites) - len(hy2.site.unique())) + " sites are not in Hydrotel")
last1 = hy2.groupby('site').last()
last_index = (last1.time == end_date)
if not all(last_index):
    print(str(sum(~last_index)) + " sites have less than a full months record")
    hy2 = hy2[hy2.site.isin(last_index.index[last_index])]
hy2 = hy2[hy2.time != end_date]
hy3 = grp_ts_agg(hy2, 'site', 'time', 'M').sum().reset_index()
hy_precip = hy3.copy()
hy_precip['mtype'] = 'precip'

### gw
#print('Getting HydroTel GWL Data:')
#
#hy_sites = mon_gw1.site.unique().copy()
#hy1 = rd_hydrotel(hy_sites, mtype='gwl_tel', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample='day', fun='avg')
#hy2 = hy1.reset_index()
#if len(hy2.site.unique()) != len(hy_sites):
#    print(str(len(hy_sites) - len(hy2.site.unique())) + " sites are not in Hydrotel")
#last1 = hy2.groupby('site').last()
#last_index = (last1.time == end_date)
#if not all(last_index):
#    print(str(sum(~last_index)) + " sites have less than a full months record")
#    hy2 = hy2[hy2.site.isin(last_index.index[last_index])]
#hy2 = hy2[hy2.time != end_date]
#hy3 = grp_ts_agg(hy2, 'site', 'time', 'M').median().reset_index()
#hy_gw = hy3.copy()
#hy_gw['mtype'] = 'gw'

### combine data and update sites

#hy_summ = concat([hy_flow, hy_precip, hy_gw]).reset_index(drop=True)
hy_summ = concat([hy_flow, hy_precip]).reset_index(drop=True)
hy_sites = hy_summ.site.unique()
mon_summ = mon_summ[mon_summ.site.isin(hy_sites)]

#################################################
#### Estimate the catchment area weights

### SW
site_catch1 = rec_catch[rec_catch.site.isin(hy_summ.site[hy_summ.mtype == 'flow'].unique())]

overlay1 = spatial_overlays(site_catch1, sw_zones, how='intersection')

overlay2 = overlay1.merge(sites1, on=['site', 'zone']).drop(['idx1', 'idx2', 'NZREACH'], axis=1)
overlay2['area'] = overlay2.area

zone_sum1 = overlay2.groupby(['zone']).area.transform('sum')
overlay2['agg_area'] = zone_sum1

overlay3 = overlay2.set_index('site').drop('geometry', axis=1)
sw_area_weight = (overlay3.area / overlay3.agg_area)
sw_area_weight.name = 'sw_area_weight'

sw_site_zone = overlay3[['zone']].copy()

## Simplify catchment polygons for eventual plotting
site_catch2 = site_catch1.drop('NZREACH', axis=1).copy()
site_catch2['geometry'] = site_catch1.simplify(100)

### precip
precip_sites1 = precip_sites[precip_sites.site.isin(hy_summ.site[hy_summ.mtype == 'precip'].unique())]
precip_site_zone = sjoin(precip_sites1, precip_zones)
precip_site_zone['precip_area_weight'] = 1/precip_site_zone.groupby(['zone'])['site'].transform('count')
precip_area_weight = precip_site_zone[['site', 'precip_area_weight']].sort_values('site').set_index('site')['precip_area_weight']

precip_site_zone1 = precip_site_zone[['site', 'zone']].set_index('site').copy()

### gw
#gw_sites1 = gw_sites[gw_sites.site.isin(hy_summ.site[hy_summ.mtype == 'gw'].unique())]
#gw_site_zone = sjoin(gw_sites1, gw_zones)
#gw_site_zone = gw_site_zone.rename(columns={'zone_left': 'zone'}).drop('zone_right', axis=1)
#gw_site_zone['gw_area_weight'] = 1/gw_site_zone.groupby(['zone'])['site'].transform('count')
#gw_area_weight = gw_site_zone[['site', 'gw_area_weight']].sort_values('site').set_index('site')['gw_area_weight']
#
#gw_site_zone1 = gw_site_zone[['site', 'zone']].set_index('site').copy()

### Combine

#area_weights = concat([sw_area_weight, precip_area_weight, gw_area_weight])
area_weights = concat([sw_area_weight, precip_area_weight])
area_weights.name = 'area_weights'

#site_zones = concat([sw_site_zone, precip_site_zone1, gw_site_zone1])
site_zones = concat([sw_site_zone, precip_site_zone1])

##############################################
#### Run the monthly stats comparisons

print('Calculating the percentiles')

def row_perc(x, mon_summ):
    mon1 = x.time.month
    mon_val = mon_summ[(mon_summ.site == x.site) & (mon_summ.mon == mon1) & (mon_summ.mtype == x.mtype)].data.values
    perc1 = percentileofscore(mon_val, x.value)
    return(perc1)

hy_summ['perc_temp'] = hy_summ.apply(row_perc, mon_summ=mon_summ, axis=1)

##############################################
#### Calc zone stats and apply categories

### Catchment area weights for SW
hy_summ1 = merge(hy_summ, area_weights.reset_index(), on='site', how='left')
if sum(hy_summ1.area_weights.isnull()) > 0:
    raise ValueError('Missing some site area weights!')
hy_summ1['perc'] = (hy_summ1['perc_temp'] * hy_summ1['area_weights']).round(2)
hy_summ1.loc[:, 'time'] = hy_summ1.loc[:, 'time'].dt.strftime('%Y-%m')

hy_summ2 = merge(hy_summ1[['mtype', 'site', 'time', 'perc']], site_zones.reset_index(), on='site', how='left')

prod1 = [sw_zones.zone.values, hy_summ2.mtype.unique(), hy_summ2.time.unique()]
mindex = MultiIndex.from_product(prod1, names=[u'zone', u'mtype', u'time'])
blank1 = Series(nan, index=mindex, name='temp')
zone_stats1 = hy_summ2.groupby(['zone', 'mtype', 'time'])['perc'].sum().round(1)
zone_stats2 = concat([zone_stats1, blank1], axis=1).perc
zone_stats2[zone_stats2.isnull()] = -1

cat_val_lst = [-10, -0.5, 10, 25, 75, 90, 100]
cat_name_lst = ['No data', 'Very low', 'Below average', 'Average', 'Above average', 'Very high']

cat1 = cut(zone_stats2, cat_val_lst, labels=cat_name_lst).astype('str')
cat1.name = 'category'
cat2 = concat([zone_stats2, cat1], axis=1)
cat3 = cat2.sort_values('perc', ascending=False).category

## Cat for catchemnts
catch_stats1 = hy_summ1.groupby(['site', 'mtype', 'time'])['perc_temp'].sum().round(1)
catch_cat1 = cut(catch_stats1, cat_val_lst, labels=cat_name_lst).astype('str')
catch_cat1.name = 'category'
catch_cat2 = concat([catch_stats1, catch_cat1], axis=1)


################################################
#### Output stats

print('Exporting results to csv')

ts_out1 = hy_summ1.copy()
ts_out2 = ts_out1.pivot_table('perc_temp', ['mtype', 'site'], 'time').round(2)
ts_out3 = ts_out2.reset_index()
ts_out3.to_csv(path.join(base_dir, ts_out_csv), index=False)

#################################################
#### Plotting

### Extract x and y data for plotting

print('Creating the plot')

def getPolyCoords(row, coord_type, geom='geometry'):
    """Returns the coordinates ('x' or 'y') of edges of a Polygon exterior"""

    # Parse the exterior of the coordinate
    exterior = row[geom].exterior

    if coord_type == 'x':
        # Get the x coordinates of the exterior
        return list(exterior.coords.xy[0])
    elif coord_type == 'y':
        # Get the y coordinates of the exterior
        return list(exterior.coords.xy[1])

zones1 = multipoly_to_poly(view_zones)

zones1['x'] = zones1.apply(getPolyCoords, coord_type='x', axis=1)
zones1['y'] = zones1.apply(getPolyCoords, coord_type='y', axis=1)

zones2 = zones1.drop('geometry', axis=1)

cant1 = GeoDataFrame(['Canterbury'], geometry=[zones1.unary_union])
cant1.columns = ['site', 'geometry']
cant1['x'] = cant1.apply(getPolyCoords, coord_type='x', axis=1)
cant1['y'] = cant1.apply(getPolyCoords, coord_type='y', axis=1)

cant2 = cant1.drop('geometry', axis=1)

## Catchments
catch1 = multipoly_to_poly(site_catch2)
catch1['x'] = catch1.apply(getPolyCoords, coord_type='x', axis=1)
catch1['y'] = catch1.apply(getPolyCoords, coord_type='y', axis=1)
catch2 = catch1.drop('geometry', axis=1)

### Combine with time series data
data1 = merge(cat1.unstack('time').reset_index(), zones2, on=['zone'])
time_index = hy_summ2.time.unique().tolist()
data1['cat'] = data1[time_index[-1]]

catch_data1 = merge(catch_cat1.unstack('time').reset_index(), catch2, on=['site'])
catch_data1['cat'] = catch_data1[time_index[-1]]

### Extract the mtype dataframes
flow_b = data1.loc[data1.mtype == 'flow']
precip_b = data1.loc[data1.mtype == 'precip']
#gw_b = data1.loc[data1.mtype == 'gw']

flow_b_catch = catch_data1.loc[catch_data1.mtype == 'flow']

flow_source = ColumnDataSource(flow_b)
precip_source = ColumnDataSource(precip_b)
#gw_source = ColumnDataSource(gw_b)
time_source = ColumnDataSource(DataFrame({'index': time_index}))

flow_catch_source = ColumnDataSource(flow_b_catch)

cant_source = ColumnDataSource(cant2)

### Set up plotting parameters
c1 = brewer['RdBu'][5]
grey1 = brewer['Greys'][7][5]

factors = cat_name_lst[::-1]
color_map = CategoricalColorMapper(factors=factors, palette=[c1[0], c1[1], c1[2], c1[3], c1[4], grey1])

### Set up dummy source for the legend
dummy_b = flow_b[['zone', 'cat', 'x', 'y']].sort_values('zone')
dummy_b.loc[:, 'cat'].iloc[0:len(factors)] = factors
dummy_source = ColumnDataSource(dummy_b)

TOOLS = "pan,wheel_zoom,reset,hover,save"

w = 700
h = w

### All three parameters

output_file(path.join(output_dir, precip_sw1_html))

## dummy figure - for legend consistency
p0 = figure(title='dummy Index', tools=[], logo=None, height=h, width=w)
p0.patches('x', 'y', source=dummy_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p0.renderers = [i for i in p0.renderers if (type(i) == renderers.GlyphRenderer) | (type(i) == annotations.Legend)]
p0.renderers[1].visible = False

## Figure 1 - precip
p1 = figure(title='Precipitation Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom', plot_height=h, plot_width=w)
p1.patches('x', 'y', source=precip_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, fill_alpha=1)
p1.renderers.extend(p0.renderers)
#p1.legend = p0.legend
p1.legend.location = 'top_left'
#Legend(p0)
hover1 = p1.select_one(HoverTool)
hover1.point_policy = "follow_mouse"
hover1.tooltips = [("Category", "@cat"), ("Zone", "@zone")]

callback1 = CustomJS(args=dict(source=precip_source), code="""
    var data = source.data;
    var f = cb_obj.value;
    source.data.cat = data[f];
    source.change.emit();
""")

select1 = Select(title='Month', value=time_index[-1], options=time_index)
select1.js_on_change('value', callback1)
#slider = Slider(start=0, end=len(time_index)-1, value=0, step=1)
#slider.js_on_change('value', callback)

layout1 = column(p1, select1)
tab1 = Panel(child=layout1, title='Precip')

## Figure 2 - flow
p2 = figure(title='Surface Water Flow Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom', plot_height=h, plot_width=w)
p2.patches('x', 'y', source=flow_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p2.renderers.extend(p0.renderers)
p2.legend.location = 'top_left'

hover2 = p2.select_one(HoverTool)
hover2.point_policy = "follow_mouse"
hover2.tooltips = [("Category", "@cat"), ("Zone", "@zone")]

callback2 = CustomJS(args=dict(source=flow_source), code="""
    var data = source.data;
    var f = cb_obj.value;
    source.data.cat = data[f];
    source.change.emit();
""")

select2 = Select(title='Month', value=time_index[-1], options=time_index)
select2.js_on_change('value', callback2)
#slider = Slider(start=0, end=len(time_index)-1, value=0, step=1)
#slider.js_on_change('value', callback)

layout2 = column(p2, select2)
tab2 = Panel(child=layout2, title='SW Flow')

## Figure 3 - GW
#p3 = figure(title='Groundwater Level Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom', plot_height=h, plot_width=w)
#p3.patches('x', 'y', source=gw_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
#p3.renderers.extend(p0.renderers)
#p3.legend.location = 'top_left'
#
#hover3 = p3.select_one(HoverTool)
#hover3.point_policy = "follow_mouse"
#hover3.tooltips = [("Category", "@cat"), ("Zone", "@zone")]
#
#callback3 = CustomJS(args=dict(source=gw_source), code="""
#    var data = source.data;
#    var f = cb_obj.value;
#    source.data.cat = data[f];
#    source.change.emit();
#""")
#
#select3 = Select(title='Month', value=time_index[-1], options=time_index)
#select3.js_on_change('value', callback3)
##slider = Slider(start=0, end=len(time_index)-1, value=0, step=1)
##slider.js_on_change('value', callback)
#
#layout3 = column(p3, select3)
#tab3 = Panel(child=layout3, title='GW Level')

## Combine
#tabs = Tabs(tabs=[tab1, tab2, tab3])
#
#show(tabs)


### Only precip and SW

#output_file(path.join(output_dir, precip_sw1_html))

tabs_alt = Tabs(tabs=[tab1, tab2])

show(tabs_alt)


### Plots for catchments

output_file(path.join(output_dir, sw1_html))

## Base figure for Canterbury


## Figure 1 - flow
p3 = figure(title='Surface Water Flow Index by catchment', tools=TOOLS, logo=None, active_scroll='wheel_zoom', plot_height=h, plot_width=w)
p3.patches('x', 'y', source=cant_source, fill_color='white', line_color="black", line_width=1)
p3.patches('x', 'y', source=flow_catch_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p3.renderers.extend(p0.renderers)
p3.legend.location = 'top_left'

hover3 = p3.select_one(HoverTool)
hover3.point_policy = "follow_mouse"
hover3.tooltips = [("Category", "@cat"), ("Site", "@site")]

callback3 = CustomJS(args=dict(source=flow_catch_source), code="""
    var data = source.data;
    var f = cb_obj.value;
    source.data.cat = data[f];
    source.change.emit();
""")

select3 = Select(title='Month', value=time_index[-1], options=time_index)
select3.js_on_change('value', callback3)
#slider = Slider(start=0, end=len(time_index)-1, value=0, step=1)
#slider.js_on_change('value', callback)

layout3 = column(p3, select3)

show(layout3)




#############################################
#### Print where results are saved

print('########################')
print('Results were saved here: ' + path.join(output_dir, ts_out_csv))
print('The plot was saved here: ' + path.join(output_dir, precip_sw1_html))




