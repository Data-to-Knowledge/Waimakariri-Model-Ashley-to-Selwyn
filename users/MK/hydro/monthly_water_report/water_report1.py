# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 16:27:09 2017

@author: MichaelEK
"""

from geopandas import read_file, overlay, sjoin, GeoDataFrame
from core.classes.hydro import hydro, all_mtypes
from core.ecan_io import rd_hydrotel
from datetime import datetime
from pandas import DateOffset, to_datetime, date_range, read_csv, concat, merge, cut, DataFrame
from os.path import join
from core.spatial.vector import spatial_overlays, multipoly_to_poly
from core.ts import grp_ts_agg, w_resample
from datetime import date
from scipy.stats import percentileofscore
from numpy import in1d, round
from bokeh.plotting import figure, save, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, LogColorMapper, Legend, CategoricalColorMapper, CustomJS
from bokeh.palettes import RdYlBu11 as palette
from bokeh.palettes import brewer
from bokeh.models.widgets import Panel, Tabs, Slider, Select
from bokeh.models.tools import WheelZoomTool
from collections import OrderedDict
from bokeh.layouts import widgetbox, column
from core.ts.met import precip_stats


###################################################
#### Parameters

base_dir = r'P:\Surface Water Quantity\Projects\Freshwater Report'
sw_poly_shp = 'sw_boundary_v01.shp'
precip_poly_shp = 'precip_boundary_v01.shp'
gw_poly_shp = 'precip_boundary_v01.shp'
rec_catch_shp = r'E:\ecan\shared\GIS_base\vector\catchments\catch_delin_recorders.shp'
streams_shp = r'E:\ecan\shared\GIS_base\vector\streams\river-environment-classification-canterbury-2010.shp'
rec_sites_shp = r'E:\ecan\shared\GIS_base\vector\catchments\recorder_sites_REC.shp'
rec_sites_details_shp = r'E:\ecan\shared\GIS_base\vector\catchments\recorder_sites_REC_details.shp'
gw_sites_shp = r'P:\Surface Water Quantity\Projects\Freshwater Report\gw_sites1.shp'

qual_codes = [10, 18, 20, 30, 50, 11, 21, 40]

month_names = ['Jan', 'Feb', 'March', 'April', 'May', 'June', 'July', 'August', 'Sept', 'Oct', 'Nov', 'Dec']

lon_zone_names = {'L': 'Lowlands', 'F': 'Foothills', 'M': 'Mountains', 'BP': 'Banks Peninsula'}

std_cat = [0.1, 1]

pot_sw_site_list_csv = 'potential_sw_site_list.csv'

### Output
catch_shp = 'recorder_catch.shp'
pot_sites_shp = 'potential_sites.shp'
pot_catch_sites_shp = 'potential_catch_sites.shp'
sw_zone_stats_shp = 'sw_zone_stats.shp'
precip_site_shp = 'precip_sites.shp'
precip_zone_stats_shp = 'precip_zone_stats.shp'

## plots
test1_html = r'E:\ecan\git\ecan_python_courses\docs\test1.html'
test2_html = r'E:\ecan\git\ecan_python_courses\docs\test2.html'

##################################################
#### Read in data

### SW
sw_poly = read_file(join(base_dir, sw_poly_shp))[['lat_zone', 'lon_zone', 'geometry']]
rec_catch = read_file(rec_catch_shp)
streams = read_file(streams_shp)
rec_sites = read_file(rec_sites_shp)
site_list = read_csv(join(base_dir, pot_sw_site_list_csv))

sw_list = site_list.replace({'lon_zone': lon_zone_names})
sw_list['zone'] = sw_list['lat_zone'] + ' - ' + sw_list['lon_zone']
sw_list = sw_list.drop(['lon_zone', 'lat_zone'], axis=1)

sw_zones = sw_poly.replace({'lon_zone': lon_zone_names})
sw_zones['zone'] = sw_zones['lat_zone'] + ' - ' + sw_zones['lon_zone']
sw_zones = sw_zones.drop(['lon_zone', 'lat_zone'], axis=1)
sw_zones['mtype'] = 'flow'

### precip
precip_sites = read_file(join(base_dir, precip_site_shp))
precip_zones = read_file(join(base_dir, precip_poly_shp))

precip_zones = precip_zones.replace({'lon_zone': lon_zone_names})
precip_zones['zone'] = precip_zones['lat_zone'] + ' - ' + precip_zones['lon_zone']
precip_zones = precip_zones.drop(['lon_zone', 'lat_zone'], axis=1)
precip_zones['mtype'] = 'precip'

### gw
gw_sites = read_file(join(base_dir, gw_sites_shp))
gw_zones = read_file(join(base_dir, gw_poly_shp))

gw_zones = gw_zones.replace({'lon_zone': lon_zone_names})
gw_zones['zone'] = gw_zones['lat_zone'] + ' - ' + gw_zones['lon_zone']
gw_zones = gw_zones.drop(['lon_zone', 'lat_zone'], axis=1)
gw_zones['mtype'] = 'gw'

gw_sites = gw_sites.replace({'lon_zone': lon_zone_names})
gw_sites['zone'] = gw_sites['lat_zone'] + ' - ' + gw_sites['lon_zone']
gw_sites = gw_sites.drop(['lon_zone', 'lat_zone'], axis=1)
#gw_sites['mtype'] = 'gw'


### Combine

zones = concat([sw_zones, precip_zones, gw_zones]).reset_index(drop=True)

#################################################
#### Select sites

### SW
sites1 = sw_list[sw_list.Notes.isnull()].drop('Notes', axis=1)

flow1 = hydro().get_data('flow', sites1.site, qual_codes)
stats_flow = flow1.stats('flow')

### precip
precip1 = hydro().get_data(mtypes='precip', sites=precip_sites.site, qual_codes=qual_codes)

### GW
gw1 = hydro().get_data(mtypes='gwl', sites=gw_sites.site, qual_codes=qual_codes)


#################################################
#### Estimate the catchment area weights

### SW
site_catch1 = rec_catch[rec_catch.site.isin(sites1.site)]

overlay1 = spatial_overlays(site_catch1, sw_zones, how='intersection')

overlay2 = overlay1.merge(sites1, on=['site', 'zone']).drop(['idx1', 'idx2', 'NZREACH'], axis=1)
overlay2['area'] = overlay2.area

zone_sum1 = overlay2.groupby(['zone']).area.transform('sum')
overlay2['agg_area'] = zone_sum1

overlay3 = overlay2.set_index('site').drop('geometry', axis=1)
sw_area_weight = (overlay3.area / overlay3.agg_area)
sw_area_weight.name = 'sw_area_weight'

sw_site_zone = overlay3[['zone']].copy()

### precip
precip_site_zone = sjoin(precip_sites, precip_zones)
precip_site_zone['precip_area_weight'] = 1/precip_site_zone.groupby(['zone'])['site'].transform('count')
precip_area_weight = precip_site_zone[['site', 'precip_area_weight']].sort_values('site').set_index('site')['precip_area_weight']

precip_site_zone1 = precip_site_zone[['site', 'zone']].set_index('site').copy()

### gw
gw_site_zone = sjoin(gw_sites, gw_zones)
gw_site_zone = gw_site_zone.rename(columns={'zone_left': 'zone'}).drop('zone_right', axis=1)
gw_site_zone['gw_area_weight'] = 1/gw_site_zone.groupby(['zone'])['site'].transform('count')
gw_area_weight = gw_site_zone[['site', 'gw_area_weight']].sort_values('site').set_index('site')['gw_area_weight']

gw_site_zone1 = gw_site_zone[['site', 'zone']].set_index('site').copy()

### Combine

area_weights = concat([sw_area_weight, precip_area_weight, gw_area_weight])
area_weights.name = 'area_weights'

site_zones = concat([sw_site_zone, precip_site_zone1, gw_site_zone1])

#################################################
#### Run monthly summary stats

### SW
flow2 = flow1.sel_ts(mtypes='flow')
flow2.index = flow2.index.droplevel('mtype')
flow3 = flow2.reset_index()

mon_flow1 = grp_ts_agg(flow3, 'site', 'time', 'M', 'median')
mon_flow1['mon'] = mon_flow1.time.dt.month
mon_flow1['mtype'] = 'flow'

### precip
precip2 = precip1.sel_ts(mtypes='precip')
precip2.index = precip2.index.droplevel('mtype')
precip3 = precip2.reset_index()

mon_precip1 = grp_ts_agg(precip3, 'site', 'time', 'M', 'sum')
mon_precip1['mon'] = mon_precip1.time.dt.month
mon_precip1['mtype'] = 'precip'

### gw
gw2 = gw1.sel_ts(mtypes='gwl')
gw2.index = gw2.index.droplevel('mtype')
gw3 = gw2.reset_index()

mon_gw1 = grp_ts_agg(gw3, 'site', 'time', 'M', 'sum')
mon_gw1['mon'] = mon_gw1.time.dt.month
mon_gw1['mtype'] = 'gw'

### Combine all mtypes

mon_summ = concat([mon_flow1, mon_precip1, mon_gw1]).reset_index(drop=True)

###############################################
#### Pull out recent monthly data from hydrotel

now1 = to_datetime(date.today())
start_date = now1 - DateOffset(months=7) - DateOffset(days=now1.day - 1)
end_date = now1 - DateOffset(days=now1.day - 1)

### SW
sites2 = sites1.copy()
sites2.loc[sites2.site.isin([64610, 65104, 68526]), 'site'] = [164610, 165104, 168526]

hy1 = rd_hydrotel(sites2.site, mtype='flow_tel', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample='day', fun='avg')
hy2 = hy1.reset_index()
if len(hy2.site.unique()) != len(sites2):
    raise ValueError("Didn't get all sites")
hy2 = hy2[hy2.time != end_date]
hy3 = grp_ts_agg(hy2, 'site', 'time', 'M', 'median')
#hy3.columns = ['site', 'mon_median_flow']

hy3.loc[:, 'site'] = hy3.site.replace([164610, 165104, 168526], [64610, 65104, 68526])

hy_flow = hy3.copy()
hy_flow['mtype'] = 'flow'

### precip
hy1 = rd_hydrotel(mon_precip1.site.unique(), mtype='precip_tel', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample='day', fun='sum')
hy2 = hy1.reset_index()
if len(hy2.site.unique()) != len(mon_precip1.site.unique()):
    raise ValueError("Didn't get all sites")
hy2 = hy2[hy2.time != end_date]
hy3 = grp_ts_agg(hy2, 'site', 'time', 'M', 'sum')
hy_precip = hy3.copy()
hy_precip['mtype'] = 'precip'

### gw
hy1 = rd_hydrotel(mon_gw1.site.unique(), mtype='gwl_tel', from_date=start_date.strftime('%Y-%m-%d'), to_date=end_date.strftime('%Y-%m-%d'), resample='day', fun='sum')
hy2 = hy1.reset_index()
if len(hy2.site.unique()) != len(mon_gw1.site.unique()):
    raise ValueError("Didn't get all sites")
hy2 = hy2[hy2.time != end_date]
hy3 = grp_ts_agg(hy2, 'site', 'time', 'M', 'sum')
hy_gw = hy3.copy()
hy_gw['mtype'] = 'gw'

### combine data

hy_summ = concat([hy_flow, hy_precip]).reset_index(drop=True)

##############################################
#### Run the monthly stats comparisons

#time_index = hy_summ.time.unique()


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
hy_summ1['perc'] = hy_summ1['perc_temp'] * hy_summ1['area_weights']

hy_summ2 = merge(hy_summ1[['mtype', 'site', 'time', 'perc']], site_zones.reset_index(), on='site', how='left')
hy_summ2.loc[:, 'time'] = hy_summ2.loc[:, 'time'].dt.strftime('%b %Y')

zone_stats1 = hy_summ2.groupby(['zone', 'mtype', 'time'])['perc'].sum().round(1)

cat_val_lst = [0, 10, 40, 60, 90, 100]
cat_name_lst = ['very low', 'below average', 'average', 'above average', 'very high']

cat1 = cut(zone_stats1, cat_val_lst, labels=cat_name_lst).astype('str')
cat1.name = 'category'
cat2 = concat([zone_stats1, cat1], axis=1)
cat3 = cat2.sort_values('perc', ascending=False).category

#################################################
#### Plotting

### Extract x and y data for plotting


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

zones1 = multipoly_to_poly(zones)

zones1['x'] = zones1.apply(getPolyCoords, coord_type='x', axis=1)
zones1['y'] = zones1.apply(getPolyCoords, coord_type='y', axis=1)

zones2 = zones1.drop('geometry', axis=1)

### Combine with time series data
data1 = merge(cat1.unstack('time').reset_index(), zones2, on=['mtype', 'zone'])
time_index = hy_summ2.time.unique().tolist()
data1['cat'] = data1[time_index[-1]]

### Extract the mtype dataframes
flow_b = data1.loc[data1.mtype == 'flow']
precip_b = data1.loc[data1.mtype == 'precip']

flow_source = ColumnDataSource(flow_b)
precip_source = ColumnDataSource(precip_b)
time_source = ColumnDataSource(DataFrame({'index': time_index}))

### Set up plotting parameters
c1 = brewer['RdBu'][5]

factors = ['very high', 'above average', 'average', 'below average', 'very low']
color_map = CategoricalColorMapper(factors=factors, palette=[c1[0], c1[1], c1[2], c1[3], c1[4]])

TOOLS = "pan,wheel_zoom,reset,hover,save"

### Plot
#output_file(test1_html)
#
#p1 = figure(title='Precipitation Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom')
#p1.patches('x', 'y', source=precip_source, fill_color={'field': 'precip_cat', 'transform': color_map}, line_color="black", line_width=1, legend='precip_cat')
#p1.legend.location = 'top_left'
##p1.toolbar.active_scroll = WheelZoomTool()
#hover1 = p1.select_one(HoverTool)
#hover1.point_policy = "follow_mouse"
#hover1.tooltips = [("Category", "@precip_cat"), ("Percentile", "@mon_precip{1.1}" + "%")]
#tab1 = Panel(child=p1, title='Precip')
#
#p2 = figure(title='Flow Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom')
#p2.patches('x', 'y', source=flow_source, fill_color={'field': 'flow_categ', 'transform': color_map}, line_color="black", line_width=1, legend='flow_categ')
#p2.legend.location = 'top_left'
##p2.toolbar.active_scroll = WheelZoomTool()
#hover2 = p2.select_one(HoverTool)
#hover2.point_policy = "follow_mouse"
#hover2.tooltips = [("Category", "@flow_categ"), ("Percentile", "@mon_flow_p{1.1}" + "%")]
#tab2 = Panel(child=p2, title='Flow')
#
#tabs = Tabs(tabs=[tab1, tab2])
#
#show(tabs)


output_file(test1_html)

## Figure 1 - precip
p1 = figure(title='Precipitation Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom')
p1.patches('x', 'y', source=precip_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p1.legend.location = 'top_left'

hover1 = p1.select_one(HoverTool)
hover1.point_policy = "follow_mouse"
hover1.tooltips = [("Category", "@cat"), ("Zone", "@zone")]

#it1 = [(i, [r1]) for i in factors]
#l1 = Legend(items=it1, location=(0, -30))
#
#p1.add_layout(l1, 'right')

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
p2 = figure(title='Flow Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom')
p2.patches('x', 'y', source=flow_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p2.legend.location = 'top_left'

hover2 = p2.select_one(HoverTool)
hover2.point_policy = "follow_mouse"
hover2.tooltips = [("Category", "@cat"), ("Zone", "@zone")]

#it1 = [(i, [r1]) for i in factors]
#l1 = Legend(items=it1, location=(0, -30))
#
#p1.add_layout(l1, 'right')

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
tab2 = Panel(child=layout2, title='Flow')

## Combine
tabs = Tabs(tabs=[tab1, tab2])

show(tabs)
























##################################################
#### Testing

gw_sites_shp = r'P:\Surface Water Quantity\Projects\Freshwater Report\gw_sites.shp'

gw1 = hydro().get_data(mtypes='gwl', sites=join(base_dir, gw_poly_shp), qual_codes=qual_codes)

gw2 = gw1.data
gw2.index = gw2.index.droplevel('mtype')
gw_stats = precip_stats(gw2)

gw_stats2 = gw_stats[(gw_stats['End time'] > '2017-01-01') & (gw_stats['Tot data yrs'] > 10)]

gw_geo1 = gw1.geo_loc.copy()

gw_sites1 = gw_stats2.index

gw_geo2 = gw_geo1.loc[gw_sites1]
gw_geo2.reset_index().to_file(gw_sites_shp)


mis_sites = mon_precip1.site.unique()[~in1d(mon_precip1.site.unique(), hy1.sites)]




t1 = grp.copy()
min1 = t1.min()
max1 = t1.max()

f1 = round((t1.values - min1)/(max1 - min1) *100, 1)

f3 = []
for i in t1.values:
    f2 = percentileofscore(t1.values, i)
    f3.append(f2)


df1 = DataFrame([t1.values, f1, f3]).T
df1.columns = ['flow_data', 'fouad', 'percentile']

df1.to_csv(join(base_dir, 'test_output_71195.csv'), index=False)


sl1 = Slider(start=0, end=len(time_index)-1, value=0, step=1)


output_file(test2_html)

p1 = figure(title='Precipitation Index', tools=TOOLS, logo=None, active_scroll='wheel_zoom')
p1.patches('x', 'y', source=precip_source, fill_color={'field': 'cat', 'transform': color_map}, line_color="black", line_width=1, legend='cat')
p1.legend.location = 'top_left'
hover1 = p1.select_one(HoverTool)
hover1.point_policy = "follow_mouse"
hover1.tooltips = [("Category", "@cat"), ("Zone", "@zone")]

callback = CustomJS(args=dict(source=precip_source, index=time_source), code="""
    var data = source.data;
    var f = cb_obj.value
    var i = index.index[f]
    cat = data[i]
    source.change.emit();
""")

slider = Slider(start=0, end=len(time_index)-1, value=0, step=1)
slider.js_on_change('value', callback)

layout = column(p1, slider)

show(p1)














