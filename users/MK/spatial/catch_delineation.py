

from core.spatial import rec_catch_del
from os.path import join
from core.ecan_io import rd_sw_rain_geo
from core.misc import select_sites

###################################
#### Parameters

base_dir = r'E:\ecan\shared\projects\climate_change\waimak_gw_model'

sites_csv = 'flow_sites.csv'
sites_shp = 'sites.shp'

#rec_sites_rec_shp = 'recorder_sites_REC.shp'

catch_del_shp = 'catch_del.shp'
#catch_del_temp_shp = r'catch_del_temp.shp'

####################################
#### Create points shapefile and go into GIS and adjust the points if needed

sites = select_sites(join(base_dir, sites_csv)).astype('int32').tolist()
sites_geo = rd_sw_rain_geo(sites).reset_index()
sites_geo.to_file(join(base_dir, sites_shp))

#### MAKE SURE YOU GO INTO GIS AND ADJUST THE SITE POINTS BEFORE THE NEXT STEP!!!

###############################
#### Catchment delination from the REC streams layers

catch1 = rec_catch_del(sites_shp=join(base_dir, sites_shp), sites_col='site', catch_output=join(base_dir, catch_del_shp))

