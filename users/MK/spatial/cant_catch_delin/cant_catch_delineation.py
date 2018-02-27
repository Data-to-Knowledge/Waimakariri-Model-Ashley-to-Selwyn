# -*- coding: utf-8 -*-
"""
Script example to delineate catchments based on site locations.
"""

from pandas import read_table, DataFrame, concat, merge, Timedelta, datetime, to_datetime, DateOffset, date_range, Timestamp, read_csv, to_numeric
from core.spatial import catch_net, pts_poly_join, flow_sites_to_shp, agg_catch, arc_catch_del, arc_spatial_join, sel_sites_poly
from geopandas import read_file
from core.ts.sw import stream_nat
from core.classes.hydro import hydro
from os.path import join
from core.ecan_io import rd_sql, rd_sw_rain_geo

from core.spatial.network import find_upstream_rec, extract_rec_catch, agg_rec_catch
from core.spatial.vector import closest_line_to_pts

###################################
#### Parameters

prod_server03 = 'SQL2012PROD03'
dw_db = 'DataWarehouse'

flow_dict = {'server': prod_server03, 'db': dw_db, 'table': 'F_HY_Flow_Data', 'site_col': 'SiteNo', 'time_col': 'DateTime', 'data_col': 'Value', 'qual_col': 'QualityCode'}

flow_stmt = 'select distinct SiteNo from F_HY_Flow_Data'

catch_shp = r'E:\ecan\shared\GIS_base\vector\catchments\river-environment-classification-watershed-canterbury-2010.shp'
streams_shp2 = r'E:\ecan\shared\GIS_base\vector\streams\rec-canterbury-2010.shp'
streams_shp = r'E:\ecan\shared\GIS_base\vector\streams\river-environment-classification-canterbury-2010.shp'
base_dir = r'P:\cant_catch_delin\recorders'

rec_sites_shp = 'rec_sites.shp'
rec_sites_rec_shp = 'recorder_sites_REC.shp'


server2 = 'SQL2012PROD03'
database2 = 'LowFlows'
sites_table = 'LowFlows.dbo.vLowFlowSite'
lowflow_sites_cols = ['RefDBaseKey']

catch_del_shp = r'catch_del.shp'
catch_del_temp_shp = r'catch_del_temp.shp'

################################
#### First define the necessary sites for the delineation
#### This can come from anywhere as long as they are int flow sites

rec_sites = rd_sql(flow_dict['server'], flow_dict['db'], stmt=flow_stmt)
rec_sites_geo = rd_sw_rain_geo(rec_sites.SiteNo.tolist())
rec_sites_geo.reset_index().to_file(join(base_dir, rec_sites_shp))


###############################
#### Catchment delination from the REC streams layers

### Load new sites layer
streams = read_file(streams_shp)
pts = read_file(join(base_dir, rec_sites_shp))

pts_seg = closest_line_to_pts(pts, streams, line_site_col='NZREACH', dis=400)
nzreach = pts_seg.copy().NZREACH.unique()
pts_seg.to_file(join(base_dir, rec_sites_rec_shp))

reaches = find_upstream_rec(nzreach)

rec_catch = extract_rec_catch(reaches)
rec_catch.to_file(join(base_dir, catch_del_temp_shp))

rec_shed = agg_rec_catch(rec_catch)
rec_shed.columns = ['NZREACH', 'geometry']
rec_shed1 = rec_shed.merge(pts_seg[['site', 'NZREACH']], on='NZREACH')

rec_shed1.to_file(join(base_dir, catch_del_shp))


##############################
#### Then aggregate the catchments to get all of the upstream area for each location








streams = read_file(streams_shp2)










