# -*- coding: utf-8 -*-
"""
Script to query the usage data from hilltop xml exports.
"""

from core.ecan_io.hilltop import iter_xml_dir, data_check_fun, proc_use_data, parse_ht_xml, all_data_fun, convert_site_names
from pandas import concat, read_csv, to_datetime, DataFrame, Series, read_hdf, merge, Grouper
from os.path import join, basename
from numpy import nan, all, array_equal, inf
from core.ecan_io import rd_sql
from core.allo_use import allo_ts_apply, allo_ts_proc


#def count_missing(df, date_col, data_col):
#    """
#    Function to count the number of missing days.
#    """
#    from pandas import Timestamp
#
#    df1 = df[[date_col, data_col]].set_index(date_col)
#    first_valid = Timestamp(df1.first_valid_index())
#    last_valid = Timestamp(df1.last_valid_index())
#    tot1 = (last_valid - first_valid).days + 1
#    count1 = df1[data_col].count()
#    mis1 = int(tot1 - count1)
#    return(mis1)


def count_missing(df, data_col):
    """
    Function to count the number of missing days.
    """
    from pandas import Timestamp, Series

    first_valid = Timestamp(df.first_valid_index())
    last_valid = Timestamp(df.last_valid_index())
    tot1 = (last_valid - first_valid).days + 1
    count1 = df[data_col].count()
    mis1 = int(tot1 - count1)
    count0 = df.loc[df[data_col] == 0, data_col].count()
    ret1 = Series([mis1, count0], index=['missing', 'zeros'])
    return(ret1)


#########################################
#### Parameters

all_use_hdf = 'E:/ecan/shared/base_data/usage/usage_daily.h5'
ht_use_csv = r'E:\ecan\shared\base_data\usage\use_daily_all_waps.csv'
wus_code = 'wus_day'
allo_csv = r'E:\ecan\shared\base_data\usage\allo_gis.csv'

wus_export = 'E:/ecan/shared/base_data/usage/wus_usage_daily.h5'

status_codes = ['Terminated - Replaced', 'Issued - Active', 'Terminated - Surrendered', 'Terminated - Cancelled', 'Terminated - Expired', 'Terminated - Lapsed', 'Issued - s124 Continuance']

allo_check_csv = r'E:\ecan\shared\base_data\usage\allo_2016_checks.csv'
allo_use_comp_csv = r'E:\ecan\shared\base_data\usage\allo_use_comp_2015-2016.csv'

########################################
#### import data

allo = read_csv(allo_csv)

all_use = read_hdf(all_use_hdf)
ht_use = read_csv(ht_use_csv)
#ht_use1 = ht_use.set_index(['wap', 'date'])['usage'].sort_index()
ht_use.columns = ['wap', 'date', 'usage']
ht_use.loc[:, 'wap'] = ht_use.loc[:, 'wap'].str.replace(',', '')

#wus = rd_sql(code='wus_day')
#wus.loc[:, 'usage'] = wus.loc[:, 'usage'].round(2)
#wus.loc[:, 'wap'] = wus.loc[:, 'wap'].str.upper()
#wus.loc[:, 'date'] = to_datetime(wus.loc[:, 'date'])
#wus1 = wus.set_index(['wap', 'date'])['usage'].sort_index()
#wus1.to_hdf(wus_export, key='wus_daily_usage', mode='w')

wus = read_hdf(wus_export).reset_index()

###############################################
#### Agg and prepare data

### Prepare allocation data
allo_wap_ann1 = allo_ts_proc(allo, start='2014-07-01', end='2016-06-30', freq='A')
allo_wap_ann2 = allo_wap_ann1.groupby(['wap', 'date'])['allo'].sum().unstack()
allo_wap_ann2.columns = ['2014-2015', '2015-2016']

### Summarise water use data by water years
wus_2015_daily = wus[(wus.date >='2014-07-01') & (wus.date <='2015-06-30')].dropna()
ht_2016_daily = ht_use[(ht_use.date >= '2015-07-01') & (ht_use.date <= '2016-06-30')].dropna()

wus_2015_ann = wus_2015_daily.groupby('wap').sum()
wus_2015_ann.columns = ['wus_2014-2015']
ht_2016_ann = ht_2016_daily.groupby('wap').sum()
ht_2016_ann.columns = ['ht_2015-2016']

### Combine annual volumes

use1 = concat([wus_2015_ann, ht_2016_ann], axis=1)
use1.index.name = 'wap'
use2 = merge(use1.reset_index(), allo_wap_ann2.reset_index(), on='wap', how='left').sort_values('2015-2016', ascending=False)

ht_wus = (use2['ht_2015-2016']/use2['wus_2014-2015']).round(3)
ht_wus.loc[ht_wus == inf] = nan
wus_allo = (use2['wus_2014-2015']/use2['2014-2015']).round(3)
wus_allo.loc[wus_allo == inf] = nan
ht_allo = (use2['ht_2015-2016']/use2['2015-2016']).round(3)
ht_allo.loc[ht_allo == inf] = nan

### Count the number of missing days and zeros
ht1 = ht_2016_daily.set_index('date')
mis1a = ht1.groupby(['wap', Grouper(level='date', freq='A-JUN')]).apply(lambda x: count_missing(x, 'usage'))

wus1 = wus_2015_daily.set_index('date')
mis1b = wus1.groupby(['wap', Grouper(level='date', freq='A-JUN')]).apply(lambda x: count_missing(x, 'usage'))

##############################################
#### Checks

### Way over allocation
## Compared to allocation
over_times = 2

way_over1 = use2[ht_allo > over_times]

## compared to previous year
way_over2 = use2[ht_wus > over_times]

### Way under
## compared to previous year
min_wus_use = 4000
under_ratio = 0.1

bool1 = use2['wus_2014-2015'] >= min_wus_use
way_under1 = use2[bool1 & (ht_wus <= under_ratio)]

## Zero usage
zero_use = use2[use2['ht_2015-2016'] == 0]

### Other
## Total WAPs in 2014-2015
sum(use2['wus_2014-2015'].notnull())

## Total WAPs in 2015-2016
sum(use2['ht_2015-2016'].notnull())

## WAPs without allocation
sum(use2['2015-2016'].isnull())

#############################################
#### Imbed flags for each into allo

allo1 = allo.copy()

name1 = 'over_' + str(over_times) + '_times_allo'
index1 = allo1.wap.isin(way_over1.wap)
allo1.loc[index1, name1] = 1

name2 = 'over_' + str(over_times) + '_times_previous_year'
index2 = allo1.wap.isin(way_over2.wap)
allo1.loc[index2, name2] = 1

name3 = 'under_' + str(int(under_ratio * 100)) + '_percent_previous_year'
index3 = allo1.wap.isin(way_under1.wap)
allo1.loc[index3, name3] = 1

name4 = 'zero_usage'
index4 = allo1.wap.isin(zero_use.wap)
allo1.loc[index4, name4] = 1

allo1.to_csv(allo_check_csv, index=False)

mis1a.columns = ['2016_missing', '2016_zeros']
mis2a = mis1a.reset_index().drop('date', axis=1)

mis1b.columns = ['2015_missing', '2015_zeros']
mis2b = mis1b.reset_index().drop('date', axis=1)

use3 = merge(use2, mis2b, on='wap', how='left')
use4 = merge(use3, mis2a, on='wap', how='left')

use4.to_csv(allo_use_comp_csv, index=False)




###############################################
#### Testing



df = ht_2016_daily.loc[ht_2016_daily.wap == 'BT26/5002', ['date', 'usage']]

d1 = ht_2016_daily.groupby('wap').apply(lambda x: count_missing(x, 'date', 'usage'))

d2 = ht_2016_daily.groupby(['wap', Grouper(key='date', freq='A-JUN')]).apply(lambda x: count_missing(x, 'usage'))

ht_2016_daily.loc[:, 'date'] = to_datetime(ht_2016_daily.loc[:, 'date'])

ht1 = ht_2016_daily.set_index('date')
d2 = ht1.groupby(['wap', Grouper(level='date', freq='A-JUN')]).apply(lambda x: count_missing(x, 'usage'))




y1 = ht1[ht1.wap == 'O32/0086']






















