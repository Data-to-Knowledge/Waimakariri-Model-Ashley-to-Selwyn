# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class
from core.classes.hydro import hydro, all_mtypes
from os import path
from core.allo_use import allo_query


#################################################
#### Load data

### Parameters
grp_by = ['date']
gw_zone = ['Eyre River']
take_type = ['Take Groundwater']

mtypes6 = 'usage'
qual_codes = [10, 18, 20, 30, 50, 11]
to_date = '2016-06-30'
base_dir = r'E:\ecan\local\Projects\requests\lee_burbery\2017-11-20'
poly = 'eyre_gw_allo_zone.shp'

export_csv = 'eyre_usage_up_to_2016-06-30.csv'
export_shp = 'eyre_usage_sites_up_to_2016-06-30.shp'

### Load in data)
allo = allo_query(grp_by=grp_by, gwaz=gw_zone, take_type=take_type, debug=True)
#h4 = hydro().get_data(mtypes=mtypes6, to_date=to_date, sites=path.join(base_dir, poly))
h4 = hydro().get_data(mtypes=mtypes6, to_date=to_date, sites=allo.wap.unique().tolist())


#############################################
#### Saving data

### Export a time series csv
h4.to_csv(path.join(base_dir, export_csv), mtypes='usage', pivot=True)

### Export a shapfile
h4.to_shp(path.join(base_dir, export_shp))












