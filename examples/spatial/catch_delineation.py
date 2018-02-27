# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 08:54:54 2017

@author: MichaelEK

Example script to delineate catchments from the REC streams and catchments.
"""

from core.spatial import rec_catch_del
from os.path import join
from core.ecan_io import rd_sw_rain_geo
from core.misc import select_sites

###################################
#### Parameters

base_dir = r'\\fileservices02\ManagedShares\Data\VirtualClimate\examples'

sites_csv = 'flow_sites.csv'
sites_shp = 'sites.shp'
catch_del_shp = 'catch_del.shp'

####################################
#### Create points shapefile and go into GIS and adjust the points if needed

sites = select_sites(join(base_dir, sites_csv)).astype('int32').tolist()
sites_geo = rd_sw_rain_geo(sites).reset_index()
sites_geo.to_file(join(base_dir, sites_shp))

#### MAKE SURE YOU GO INTO GIS AND ADJUST THE SITE POINTS BEFORE THE NEXT STEP!!!

###############################
#### Catchment delination from the REC streams layers

catch1 = rec_catch_del(sites_shp=join(base_dir, sites_shp), sites_col='site', catch_output=join(base_dir, catch_del_shp))

