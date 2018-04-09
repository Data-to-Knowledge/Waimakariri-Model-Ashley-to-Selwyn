# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/04/2018 4:42 PM
"""

from __future__ import division
from core import env
from pdsql.mssql import rd_sql
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.nitrate_at_key_receptors import \
    get_well_ids, get_str_ids
import pandas as pd
import geopandas as gpd

shp_data = gpd.read_file(r"C:\Users\MattH\Downloads\all_wells_2.shp")
groups = set(shp_data.loc[shp_data.wdc_group.notnull(),'wdc_group']) - {''}

outdata = {}
for zone in groups:
    sites = shp_data.loc[shp_data.wdc_group==zone, 'Field1']

    data = rd_sql(server='sql2012dev01',
                  database='hydro',
                  table='HilltopWQMtypes',
                  col_names=['SiteID', 'MeasurementType', 'CollectionTime', 'Value'],
                  where_col={'SiteID': sites, 'MeasurementType': ['Nitrate Nitrogen', 'Nitrate-N + Nitrite-N'],
                             'Param': ['data']},
                  from_date='2014-12-31',
                  to_date='2017-12-31',
                  date_col='CollectionTime')
    idx = data.Value.str.contains('<')
    data.loc[:, 'Value'] = data.Value.str.replace('<', '').astype(float)
    data.loc[idx, 'Value'] *= 0.5
    outdata[zone] = data.describe()['Value']
outdata = pd.DataFrame(outdata).transpose()
outdata.to_csv(r"T:\Temp\temp_gw_files\missing_wdc_n_data_pull.csv")
