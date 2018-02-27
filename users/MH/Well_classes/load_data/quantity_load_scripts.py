"""
Author: matth
Date Created: 28/02/2017 11:50 AM
"""

from __future__ import division
import pandas as pd
from core.ecan_io import rd_sql, sql_db

daily_quant_type_codes = {
    'A': 'Aquifer Test',
    'D': 'Daily Instantaneous Recorded (noon)',
    'F': 'Owner - Not Quality Assured',
    'N': 'Manual',
    'O': 'Owner - Quality Assured',
    'X': 'Suspect Data',
    'Y': 'Weekly Instantaneous Recorded (Wed)'
}

# debug
# the following load scripts bring in the well level data to the well quant item
def load_yearly(well_list, df=None, start_date=None, end_date=None):
    if df is None:
        df = load_yearly_quant_df(well_list)
    df = df[df.WELL_NO.isin(well_list)]
    df = df.drop(['WELL_NO'], axis=1)

    if start_date is not None and end_date is not None:
        df = df[(df.quant_time >= start_date) & (df.quant_time <= end_date)]

    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_yearly_quant_df(well_list):
    data_names = {'Year_Read': 'quant_year',
                  'AVERAGE_WATER_LEVEL': 'water_level',
                  'Count': 'num_samp',
                  'Max_Water_Level': 'max_water_level',
                  'Min_Water_Level': 'min_water_level'}
    alldata = rd_sql(**sql_db.wells_db.yearly)
    df = alldata[alldata.WELL_NO.isin(well_list)]

    df.rename(columns=data_names)
    temp_time = ['01/06/{}'.format(e) for e in df.quant_year]
    df['quant_time'] = temp_time
    df.quant_time = pd.to_datetime(df.quant_time, format='%d/%m/%Y')
    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_monthly(well_list, df=None, start_date=None, end_date=None):
    """
    :param well_list:
    :param df:
    :param start_date:
    :param end_date:
    :return:
    """
    if df is None:
        df = load_monthly_quant_df(well_list)
    df = df[df.WELL_NO.isin(well_list)]
    df = df.drop(['WELL_NO'], axis=1)

    if start_date is not None and end_date is not None:
        df = df[(df.quant_time >= start_date) & (df.quant_time <= end_date)]

    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_monthly_quant_df(well_list):
    data_names = {'DMY': 'quant_time',
                  'Year': 'quant_year',
                  'Month': 'quant_month',
                  'AverageWL': 'water_level',
                  'Count': 'num_samp',
                  'FirstDateInMonth': 'quant_first_dt',
                  'LastDateinMonth': 'quant_last_dt'}
    alldata = rd_sql(**sql_db.wells_db.monthly)
    df = alldata[alldata.WELL_NO.isin(well_list)]

    df.rename(columns=data_names)
    df.quant_time = pd.to_datetime(df.quant_time, format='%d/%m/%Y')
    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_daily(well_list, df=None, start_date=None, end_date=None):
    if df is None:
        df = load_daily_quant_df(well_list)
    df = df[df.WELL_NO.isin(well_list)]
    df = df.drop(['WELL_NO', 'DateEntered', 'LogonIdEntered', 'timestamp'], axis=1)

    if start_date is not None and end_date is not None:
        df = df[(df.quant_time >= start_date) & (df.quant_time <= end_date)]
    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_daily_quant_df(well_list):
    data_names = {'Id': 'quant_id',
                  'DEPTH_TO_WATER': 'water_level',
                  'DATE_READ': 'quant_time',
                  'TIDEDA_FLAG': 'quant_type',
                  'COMMENTS': 'comments'}
    alldata = rd_sql(**sql_db.wells_db.daily)
    df = alldata[alldata.WELL_NO.isin(well_list)]

    df.rename(columns=data_names)

    df.quant_time = pd.to_datetime(df.quant_time)
    df.quant_time[df.quant_time.dt.hour == 0] = df.quant_time[df.quant_time.dt.hour == 0] + pd.DateOffset(hours=12)

    df.sort(columns='quant_time')
    df = df.reindex(index=range(len(df.quant_time)))

    return df


def load_telemetered(df=None):
    if df is None:
        df = load_telemetered_quant_df()
    raise ValueError('not implemented')


def load_telemetered_quant_df():
    raise ValueError('not implemented')


def load_best_frequency(df=None):
    raise ValueError('not implemented')
