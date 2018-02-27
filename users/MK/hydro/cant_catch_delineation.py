# -*- coding: utf-8 -*-
"""
Script example to delineate catchments based on site locations.
"""

from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.spatial import catch_net, pts_poly_join, flow_sites_to_shp, agg_catch, arc_catch_del, arc_spatial_join
from geopandas import read_file
from core.ts.sw import stream_nat
from core.classes.hydro import hydro
from os.path import join

###################################
#### Parameters

sites_r = [69505, 69514, 69508]
sites_m = [69515, 69520, 169535]
qual_codes = [10, 18, 20, 30, 50]

streams_shp = r'S:\Surface Water\shared\GIS_base\vector\streams\rec_mfe_cant_no_1st.shp'

base_dir = r'P:\cant_catch_delin\set1'
bound_shp = 'orari_catch.shp'

catch_del_shp = r'catch_del.shp'
catch_sites_csv = r'results\catch_sites.csv'
export_sites_shp = 'sites.shp'
export_catch_shp = r'results\catch_del_poly.shp'

################################
#### First define the necessary sites for the delineation
#### This can come from anywhere as long as they are int flow sites

#### Load in data and create shapefile

h1 = hydro().get_data('flow', sites_r, qual_codes)
h2 = h1.get_data('flow_m', sites_m)

sites_geo = flow_sites_to_shp(h2.sites, export=True, export_path=join(base_dir, export_sites_shp))

###############################
#### The run the arcgis catchment delineation script in arcgis
#### Use the export from sites_geo as the sites_in path and use a polygon shpaefile that defines the desired boundary
#### Define the output working directory and run the script through!

#### Make sure to check that the sites are in the correct location!!!
#### check that after the processing that the delineations make sense!!!

arc_catch_del(base_dir, join(base_dir, bound_shp), join(base_dir, export_sites_shp), site_num_col='site', point_dis=1000, stream_depth=10, grid_size=8, pour_dis=20, streams=streams_shp, overwrite_rasters=False)
arc_spatial_join(base_dir)

##############################
#### Then aggregate the catchments to get all of the upstream area for each location

catch_del = agg_catch(join(base_dir, catch_del_shp), join(base_dir, catch_sites_csv))
catch_del.to_file(join(base_dir, export_catch_shp))

#############################
#### The result of the last function can be used to extract many types of point data
#### For example to return all abstractions within each catchment

#allo_loc = read_file(allo_loc_shp)[[crc_col, 'geometry']]
#crc_catch2, catch4 = pts_poly_join(allo_loc, catch_del, 'site', pts_id_col=crc_col)
















