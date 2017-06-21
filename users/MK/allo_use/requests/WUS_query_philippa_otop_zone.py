# -*- coding: utf-8 -*-
"""
Created on Wed Jun 08 14:39:29 2016

@author: MichaelEK
"""

from pandas import merge, read_csv, DataFrame
from query_use_allo_fun import w_query

#################################
### Parameters

series_csv = 'C:/ecan/base_data/usage/usage_takes_series_with_cav.csv'
allo_csv = 'C:/ecan/base_data/usage/takes_results_with_cav.csv'

allo_cols = ['crc', 'wap', 'take_type', 'catchment', 'irr_area', 'gw_zone', 'sw_zone', 'use_type', 'catchment_num', 'cwms_zone']
cwms_zone = ['Orari-Temuka-Opihi-Pareora']
take_type = ['Take Groundwater']

path = 'C:/ecan/Projects/requests/philippa/'

export_csv1 = 'old_gw_zones_otop_v05.csv'
export_csv2 = 'new_gw_zones_otop_v05.csv'

#export_zone_csv = 'otop_zone_tot_v04.csv'
#export_zone_type_csv = 'otop_zone_use_type_v02.csv'

#################################
### Read in allocation and usage data and merge data

series = read_csv(series_csv)
allo = read_csv(allo_csv)[allo_cols]

allo_use1 = merge(series, allo, on=['crc', 'wap'])

### Read in input data to be used in the query

sites = read_csv(path + 'input_data.csv')

### Merge new gw zones to base data

allo_use2 = merge(allo_use1, sites, on='crc', how='left')

#################################
### Query data

q1 = w_query(allo_use2, grp_by=['dates', 'gw_zone'], take_type=take_type, cwms_zone=cwms_zone, export_path=path + export_csv1)

q2 = w_query(allo_use2, grp_by=['dates', 'new_gw_zone'], take_type=take_type, cwms_zone=cwms_zone, export_path=path + export_csv2)

#q3 = w_query(allo_use2, grp_by=['dates'], take_type=take_type, cwms_zone=cwms_zone, export_path=path + export_zone_csv)
#
#q4 = w_query(allo_use2, grp_by=['dates', 'use_type'], years=[2015], take_type=take_type, cwms_zone=cwms_zone, export_path=path + export_zone_type_csv)


q10 = w_query(allo_use2, grp_by=['dates'], years=[2015], take_type=take_type, cwms_zone=cwms_zone, export_path=path + export_zone_csv, export=False, debug=True)


q10[(q10.use_type == 'public_supply') & (q10.usage_m3.notnull())]

q10[q10.crc == 'CRC101875']

allo1[allo1.crc == 'CRC101875']

takes5[takes5.crc == 'CRC101875']



ACT072759 - GW
ACT073220 - SW

allo_use2[allo_use2.crc == 'CRC101875']

'CRC052684'


