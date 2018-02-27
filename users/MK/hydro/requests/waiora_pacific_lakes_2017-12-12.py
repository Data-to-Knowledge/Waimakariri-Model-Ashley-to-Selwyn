# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 09:03:03 2017

@author: MichaelEK
"""

from core.ecan_io.hilltop import rd_squalarc
from core.classes.hydro import hydro
from core.classes.hydro.misc import mtype_codes
from core.classes.hydro.ecan_import import rd_sw_rain_geo
from core.ecan_io import rd_sql
from core.spatial.vector import sel_sites_poly, xy_to_gpd
from core.plotting.bokeh.tools import getCoords
from os import path
from pandas import to_numeric, merge, concat

###################################
### Parameters

server = 'SQL2012PROD03'
db = 'Hydstra'
period_tab = 'PERIOD'
var_tab = 'variable'

period_cols = ['STATION', 'VARFROM', 'VARIABLE', 'PERSTART', 'PEREND']
period_names = ['site', 'varfrom', 'varto', 'start', 'end']

var_cols = ['VARNUM', 'VARNAM']
var_names = ['varto', 'parameter']

base_dir = r'E:\ecan\local\Projects\requests\waiora_pacific'

poly_shp = 'lakes.shp'

mtype1 = 'lake_wl_rec_qc'

## Output
quan_csv = 'quan_sites_summary_2017-12-12.csv'
qual_csv = 'qual_sites_summary_2017-12-12.csv'

qual_data_csv = 'qual_data_2017-12-12.csv'

##################################
### Read in data

## Hydstra codes to mtypes
mtype_hydstra_codes = {1: [100], 2: [140, 143], 3: [450], 4: [110], 5: [10], 6: [130]}
#    unit_codes = {100: 'meter', 10: 'mm', 110: 'm', 140: 'm**3/s', 143: 'l/s', 450: 'degC'}
device_data_type = {100: 'mean', 140: 'mean', 143: 'mean', 450: 'mean', 110: 'mean', 130: 'mean', 10: 'tot'}

## Removals
rem_dict = {'165131': [140, 140], '69302': [140, 140], '71106': [140, 140]}

### Import
period1 = rd_sql(server, db, period_tab, period_cols, where_col='DATASOURCE', where_val=['A'])
period1.columns = period_names
period1.loc[:, 'site'] = period1.site.str.strip()

### Determine the variables to extract
period2 = period1[period1.varto.isin(period1.varto.round())].sort_values('site')
period2 = period2[period2.varto != 101]
for i in rem_dict:
    period2 = period2[~((period2.site == i) & (period2.varfrom == rem_dict[i][0]) & (period2.varto == rem_dict[i][1]))]

### Extract the sites for the specific mtypes
mtype_num = mtype_codes[mtype1]
varto_list = mtype_hydstra_codes[mtype_num]

period3 = period2[period2.varto.isin(varto_list)].reset_index(drop=True)


### GIS
sites1 = to_numeric(period2.site, 'coerce').dropna().astype('int32').tolist()

site_geo = rd_sql('SQL2012PROD05', 'Bgauging', 'RSITES', col_names=['SiteNumber', 'River', 'SiteName', 'NZTMX', 'NZTMY'], where_col='SiteNumber', where_val=sites1)

site_geo.columns = ['site', 'feature', 'name', 'NZTMX', 'NZTMY']
site_geo.loc[:, 'site'] = to_numeric(site_geo.loc[:, 'site'], errors='ignore')

site_geo2 = xy_to_gpd(df=site_geo, id_col=['site', 'feature', 'name'], x_col='NZTMX', y_col='NZTMY')
site_geo3 = site_geo2.loc[site_geo2.site > 0, :]
site_geo3.loc[:, 'site'] = site_geo3.loc[:, 'site'].astype('int32')

sites2 = sel_sites_poly(site_geo3, path.join(base_dir, poly_shp))

### Get parameters
## Quantity
param1 = period2[period2.site.isin(sites2.site.astype(str))]
#param1.sort_values('varto')

var1 = rd_sql(server, db, var_tab, var_cols, {'VARNUM': param1.varto.astype('int32').tolist()}, rename_cols=var_names)
var1['parameter'] = var1['parameter'].str.strip()

param2 = merge(param1, var1, on='varto')
param2['site'] = param2['site'].astype('int32')

sites3 = getCoords(sites2).drop('geometry', axis=1)
sites3.columns = ['site', 'feature', 'name', 'NZTMX', 'NZTMY']

param3 = merge(param2, sites3, on='site').drop(['varfrom', 'varto'], axis=1)

## Quality
wq_data1 = rd_squalarc(path.join(base_dir, poly_shp))
wq_data2 = wq_data1[wq_data1.site.str.startswith('SQ')].copy()

sw_sites_tab = rd_sql('SQL2012PROD05', 'Squalarc', 'SITES', col_names=['SITE_ID', 'SOURCE', 'SITE_NAME', 'NZTMX', 'NZTMY'])
sw_sites_tab.columns = ['site', 'feature', 'name', 'NZTMX', 'NZTMY']

xy_sw = sw_sites_tab[sw_sites_tab.site.isin(wq_data2.site)]

grp1 = wq_data2.groupby(['site', 'parameter'])

min_date = grp1['date'].min()
min_date.name = 'start'
max_date = grp1['date'].max()
max_date.name = 'end'
count1 = grp1['date'].count()
count1.name = 'sample_count'

wq_data3 = concat([min_date, max_date, count1], axis=1)

wq_data4 = merge(wq_data3.reset_index(), xy_sw, on='site')

### Export summaries
param3.to_csv(path.join(base_dir, quan_csv), index=False)
wq_data4.to_csv(path.join(base_dir, qual_csv), index=False)


wq_data2.to_csv(path.join(base_dir, qual_data_csv), index=False, encoding='utf-8')





