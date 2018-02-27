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
min_flow_fields = ['SiteID', 'BandNo', 'PeriodNo', 'isTrigger', 'Allocation', 'Flow']
site_band_fields = ['SiteID', 'BandNo', 'isActive']

min_flow_table = 'LowFlowSiteBandPeriodAllocation'
site_band_table = 'LowFlowSiteBand'

output_csv1 = r'E:\ecan\local\Projects\requests\low_flows\2017-12-07\missing_triggers.csv'
output_csv2 = r'E:\ecan\local\Projects\requests\low_flows\2017-12-07\max_allo_below_100.csv'

###########################################
### Pull in data and query

restr_val = rd_sql(server1, database1, min_flow_table, min_flow_fields)
site_band = rd_sql(server1, database1, site_band_table, site_band_fields)

grp1 = restr_val.groupby(['SiteID', 'BandNo', 'PeriodNo'])['isTrigger'].sum()
miss_sites = grp1[grp1 == 0]
miss_sites.name = 'isTriggerCount'

miss_sites1 = miss_sites.reset_index()

miss_sites2 = merge(miss_sites1, site_band, on=['SiteID', 'BandNo'], how='left')

miss_sites2.to_csv(output_csv1, header=True, index=False)


grp2 = restr_val.groupby(['SiteID', 'BandNo', 'PeriodNo'])
min1 = grp2['Allocation'].min()
max1 = grp2['Allocation'].max()

#min1[min1 >= 100]
max2 = max1[max1 < 100].copy().reset_index()
max3 = merge(max2, site_band, on=['SiteID', 'BandNo'], how='left')
max3.to_csv(output_csv2, header=True, index=False)


