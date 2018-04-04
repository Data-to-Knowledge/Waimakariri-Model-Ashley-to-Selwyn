# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/04/2018 4:42 PM
"""

from __future__ import division
from core import env
from pdsql.mssql import rd_sql

sites = ['SQ32994',
'SQ30267',
'SQ30261',
'SQ30258',]
data = rd_sql(server='sql2012dev01',
              database='hydro',
              table='HilltopWQMtypes',
              col_names=['SiteID', 'MeasurementType', 'CollectionTime', 'Value'],
              where_col={'SiteID': sites, 'MeasurementType': ['Nitrate Nitrogen', 'Nitrate-N + Nitrite-N'], 'Param': ['data']},

              to_date='2017-12-31',
              date_col='CollectionTime')
idx = data.Value.str.contains('<')
data.loc[:, 'Value'] = data.Value.str.replace('<', '').astype(float)
data.loc[idx, 'Value'] *= 0.5
data.to_csv(r"T:\Temp\temp_gw_files\n_data_pull.csv")
