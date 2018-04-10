# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/04/2018 9:35 AM
"""

from __future__ import division
from core import env
from pdsql.mssql import rd_sql
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import get_all_well_row_col
import geopandas as gpd

if __name__ == '__main__':

    all_wells = get_all_well_row_col()
    groups = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\support_shp\cl_o18_data.shp")
    sites = list(groups.loc[groups.well_group != '','Field1'])
    mapper = groups.loc[groups.well_group != ''].set_index('Field1').loc[:,'well_group'].to_dict()
    data = rd_sql(server='sql2012dev01',
                  database='hydro',
                  table='HilltopWQMtypes',
                  col_names=['SiteID', 'MeasurementType', 'Value'],
                  where_col={'SiteID':sites, 'MeasurementType': ['Oxygen 18', 'Chloride'],'Param':['data']},
                  from_date='2008-01-01',
                  to_date='2017-12-31',
                  date_col='CollectionTime',
                  ).reset_index()
    data.loc[:,'grouper'] = data.loc[:,'SiteID']
    data = data.replace({'grouper': mapper})
    idx = data.Value.str.contains('<')
    data.loc[:,'Value'] = data.Value.str.replace('<','').astype(float)
    data.loc[idx,'Value'] *= 0.5
    summary_o18 = data.loc[data.MeasurementType=='Oxygen 18',['grouper', 'Value']].groupby('grouper').describe().loc[:,'Value']
    summary_cl = data.loc[data.MeasurementType=='Chloride',['grouper', 'Value']].groupby('grouper').describe().loc[:,'Value']
    outdata = pd.merge(summary_o18, summary_cl, left_index=True, right_index=True, suffixes=('_o18', '_cl'))
    outdata.to_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\targets_for_well_adjust\cl_o18_data_grouped.csv")

