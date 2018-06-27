# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/06/2018 9:04 AM
"""

from __future__ import division
from core import env
from core.ecan_io import rd_sql, sql_db
import numpy as np
import pandas as pd
from core.classes.hydro import hydro
import rasterio
from copy import deepcopy

if __name__ == '__main__':
    elv_sheet = pd.read_excel(
        env.sci('Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/xyz.xlsx'),
        index_Col=0)
    elv_sheet = elv_sheet.set_index('well')

    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4)]  # keep only waimak zones
    well_details = well_details[pd.notnull(well_details['DEPTH'])]

    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:, 'WELL_NO'] = [e.strip() for e in screen_details.loc[:, 'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    well_details.loc[:, 'WELL_NO'] = [e.strip() for e in well_details.loc[:, 'WELL_NO']]
    well_details = well_details.set_index('WELL_NO')

    data = hydro().get_data(mtypes=['gwl'], sites=list(well_details.index), from_date='2013-01-01').data
    data = data.loc['aq_wl_cont_qc'].reset_index()
    data1 = hydro().get_data(mtypes=['gwl_m'], sites=list(well_details.index), from_date='2013-01-01').data
    temp = data1.groupby(level=['mtype', 'site']).describe()[
        ['min', '25%', '50%', '75%', 'mean', 'max', 'count']].round(2)
    sites = list(temp.loc['aq_wl_disc_qc'].index[temp.loc['aq_wl_disc_qc']['count'] >= 10])

    data1 = data1.loc['aq_wl_disc_qc', sites].reset_index()
    data = pd.concat((data, data1), axis=0)
    data['month'] = [e.month for e in data.time]
    data['year'] = [e.year for e in data.time]
    data.loc[data.data > 990, 'data'] = np.nan
    data.loc[data.data < -990, 'data'] = np.nan
    temp = data.groupby(['site', 'year'])
    out_data = pd.DataFrame(temp.max()['data'])
    out_data['count'] = temp.count()['data']
    out_data = out_data.reset_index()
    out_data.loc[:, 'nztmx'] = out_data.site.replace(well_details['NZTMX'].to_dict())
    out_data.loc[:, 'nztmy'] = out_data.site.replace(well_details['NZTMY'].to_dict())
    out_data.loc[:, 'depth'] = out_data.site.replace(well_details['DEPTH'].to_dict())
    out_data.loc[:, 'REFERENCE_RL'.lower()] = out_data.site.replace(well_details['REFERENCE_RL'].to_dict())
    out_data.loc[:, 'GROUND_RL'.lower()] = out_data.site.replace(well_details['GROUND_RL'].to_dict())
    out_data.loc[:, 'd2water'] = out_data.loc[:, 'data'] + -out_data.loc[:, 'ground_rl']
    out_data.loc[:, 'water_elv'] = out_data.loc[:, 'data'] + out_data.loc[:, 'REFERENCE_RL'.lower()]
    out_data2 = out_data.dropna()
    out_data2.drop(['data', 'reference_rl', 'ground_rl'], axis=1, inplace=True)
    out_data3 = out_data2.groupby('site').aggregate({'nztmx': 'first',
                                                     'nztmy': 'first',
                                                     'depth': 'first',
                                                     'd2water': np.min,
                                                     'water_elv': np.max})
    out_data2.to_csv(r"C:\Users\MattH\Downloads\wdc_data_pull\water_level_data.csv")
    out_data3.to_csv(r"C:\Users\MattH\Downloads\wdc_data_pull\highest_measured.csv")

