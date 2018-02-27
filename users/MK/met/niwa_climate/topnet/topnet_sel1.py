# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:38:19 2017

@author: MichaelEK
"""

from xarray import open_dataset, open_mfdataset, concat
from os import path
from core.misc import rd_dir, unarchive_dir, save_df, get_subdir
from geopandas import read_file
from core.ecan_io.met import nc_add_gis
from pandas import to_datetime, DataFrame
from core.ts.met.topnet import proc_topnet_nc
from core.spatial.vector import closest_line_to_pts
from core.ts import grp_ts_agg

###########################################
#### Parameters

base_path = r'E:\ecan\shared\base_data\niwa\climate_projections\topnet\waimak2'

rec_shp = r'E:\ecan\shared\GIS_base\vector\streams\rec-canterbury-2010.shp'

new_rec = r'E:\ecan\shared\base_data\niwa\climate_projections\topnet\waimak2\rec_reaches.shp'

param = 'flow_rate'

out_path = r'E:\ecan\shared\projects\climate_change\waimak_gw_model'
sites_shp = 'sites.shp'
rcppast_hdf = 'rcppast.h5'
rcpproj_hdf = 'rcpproj.h5'
sites_rec_shp = 'sites_rec.shp'

##########################################
#### Read data

dirs = get_subdir(base_path)
nc_file1 = path.join(base_path, dirs[0], rd_dir(path.join(base_path, dirs[0]), 'nc')[0])

nc1 = open_mfdataset(nc_file1)
rec = read_file(rec_shp)
sites = read_file(path.join(out_path, sites_shp))

#########################################
#### Select only the REC reaches with data and export new REC layer

reach1 = nc1.nrch.values
nc1.close()

rec2 = rec[rec.NZREACH.isin(reach1)]
rec2.to_file(new_rec)

#########################################
#### Find the closest reach to the points

sites_reach1 = closest_line_to_pts(sites, rec2, 'NZREACH', 400)

########################################
#### Extract reach data from topnet results

flow_dict = {}
for i in dirs:
    nc_files = rd_dir(path.join(base_path, i), 'nc')
    nc_file_paths = [path.join(base_path, i, h) for h in nc_files]
    df_temp = DataFrame()
    for nc in nc_file_paths:
        model = nc.split(i + '_')[1].split('.nc')[0]
        nc2 = open_mfdataset(nc)
        df1 = nc2.sel(nrch=sites_reach1.NZREACH.values)[param].to_dataframe()[param]
        df1.name = model
        df_temp = concat([df_temp, df1], axis=1)
    flow_dict.update({i: df_temp})

rcppast = flow_dict.pop('RCPpast')
rcpproj = concat(flow_dict, axis=1)

#######################################
#### Save results

save_df(rcppast, path.join(out_path, rcppast_hdf))
save_df(rcpproj, path.join(out_path, rcpproj_hdf))
sites_reach1.to_file(path.join(out_path, sites_rec_shp))














##############################################
#### Testing

top2 = top1.reset_index()
t1 = top2[top2.nrch == 13041967].drop('nrch', axis=1).set_index('time')


proj1 = rcpproj.copy().reset_index()

proj2 = grp_ts_agg(proj1, 'nrch', 'time', 'A-JUN').mean()






























