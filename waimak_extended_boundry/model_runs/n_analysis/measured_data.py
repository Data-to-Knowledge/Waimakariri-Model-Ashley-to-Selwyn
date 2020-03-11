# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/04/2018 10:22 AM
This script no longer works due to changes in ecan data structure, it is left here as is where is

"""

from __future__ import division
from pdsql.mssql import rd_sql
import pandas as pd


if __name__ == '__main__':
    # waiamak wells # from 2014 to present
    sites = list(pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\WaimakPlainsActiveWells.xlsx').iloc[:,0])
    data = rd_sql(server='sql2012dev01',
                  database='hydro',
                  table='HilltopWQMtypes',
                  col_names=['SiteID', 'MeasurementType', 'CollectionTime', 'Value'],
                  where_col={'SiteID':sites, 'MeasurementType':['Nitrate Nitrogen'],'Param':['data']},
                  from_date='2014-01-01',
                  to_date='2017-12-31',
                  date_col='CollectionTime')
    idx = data.Value.str.contains('<')
    data.loc[:,'Value'] = data.Value.str.replace('<','').astype(float)
    data.loc[idx,'Value'] *= 0.5
    summary = data.loc[:,['SiteID','Value']].groupby('SiteID').describe()
    data.to_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\waimak_n_2014_2017.csv")
    summary.to_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\waimak_n_2014_2017_summary.csv")

    # all canty wells # for years 2014 to end 2017
    sites = list(pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\AllActiveCantyPlainsWells.xlsx').iloc[:,0])
    data = rd_sql(server='sql2012dev01',
                  database='hydro',
                  table='HilltopWQMtypes',
                  col_names=['SiteID', 'MeasurementType', 'CollectionTime', 'Value'],
                  where_col={'SiteID':sites, 'MeasurementType':['Nitrate Nitrogen'],'Param':['data']},
                  from_date='2014-01-01',
                  to_date='2017-12-31',
                  date_col='CollectionTime')
    idx = data.Value.str.contains('<')
    data.loc[:,'Value'] = pd.to_numeric(data.Value.str.replace('<',''), errors='coerce')
    data.loc[idx,'Value'] *= 0.5
    data.loc[:,'year'] = [e.year for e in data.CollectionTime]
    summary = data.loc[:,['SiteID','Value','year']].groupby(['SiteID','year']).describe()
    data.to_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\all_canterbury_n_2014_2017.csv")
    summary.to_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\Spreadsheets\all_canterbury_n_2014_2017_summary.csv")


