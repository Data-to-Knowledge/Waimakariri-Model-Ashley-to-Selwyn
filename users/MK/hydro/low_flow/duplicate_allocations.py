# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 19:05:29 2017

@author: MichaelEK
"""

from pandas import merge, to_datetime, read_csv, concat
from core.misc.misc import printf, select_sites
from core.ecan_io import rd_sql
from core.spatial.vector import sel_sites_poly
from numpy import ndarray, in1d, nan, array
from datetime import date

############################################
### Parameters

server1 = 'SQL2012PROD03'
database1 = 'LowFlows'
min_flow_fields = ['SiteID', 'BandNo', 'PeriodNo', 'Allocation', 'Flow']
site_band_fields = ['SiteID', 'BandNo', 'isActive']

min_flow_table = 'LowFlowSiteBandPeriodAllocation'
site_band_table = 'LowFlowSiteBand'

output_csv1 = r'E:\ecan\local\Projects\requests\low_flows\2017-12-07\missing_triggers.csv'
output_csv2 = r'E:\ecan\local\Projects\requests\low_flows\2017-12-07\max_allo_below_100.csv'

###########################################
### Pull in data and query

restr_val = rd_sql(server1, database1, min_flow_table, min_flow_fields)
site_band = rd_sql(server1, database1, site_band_table, site_band_fields)

restr_val1 = restr_val.groupby(['SiteID', 'BandNo', 'PeriodNo', 'Allocation'])['Flow'].count()
restr_val2 = restr_val1[restr_val1 > 1].copy()
restr_val2.name = 'dup_allo'

bad_sites1 = merge(restr_val2.reset_index(), site_band, on=['SiteID', 'BandNo'], how='left')

miss_sites2.to_csv(output_csv1, header=True, index=False)


grp2 = restr_val.groupby(['SiteID', 'BandNo', 'PeriodNo'])
min1 = grp2['Allocation'].min()
max1 = grp2['Allocation'].max()

#min1[min1 >= 100]
max2 = max1[max1 < 100].copy().reset_index()
max3 = merge(max2, site_band, on=['SiteID', 'BandNo'], how='left')
max3.to_csv(output_csv2, header=True, index=False)


