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

if __name__ == '__main__':

    # todo make a list of sites for all wells
    all_wells = get_all_well_row_col()
    sites = all_wells.index
    data = rd_sql(server='sql2012dev01',
                  database='hydro',
                  table='HilltopWQMtypes',
                  col_names=['SiteID', 'MeasurementType', 'Value'],
                  where_col={'SiteID':sites, 'MeasurementType': ['Oxygen 18', 'Chloride'],'Param':['data']},
                  from_date='2008-01-01',
                  to_date='2017-12-31',
                  date_col='CollectionTime',
                  ).reset_index()
    idx = data.Value.str.contains('<')
    data.loc[:,'Value'] = data.Value.str.replace('<','').astype(float)
    data.loc[idx,'Value'] *= 0.5
    summary_o18 = data.loc[data.MeasurementType=='Oxygen 18',['SiteID', 'Value']].groupby('SiteID').describe().loc[:,'Value']
    summary_cl = data.loc[data.MeasurementType=='Chloride',['SiteID', 'Value']].groupby('SiteID').describe().loc[:,'Value']
    outdata = pd.merge(summary_o18, summary_cl, left_index=True, right_index=True, suffixes=('_o18', '_cl'))
    outdata = pd.merge(outdata, all_wells, left_index=True, right_index=True)
    outdata.to_csv(r"C:\Users\MattH\Downloads\cl_o18_data.csv")

