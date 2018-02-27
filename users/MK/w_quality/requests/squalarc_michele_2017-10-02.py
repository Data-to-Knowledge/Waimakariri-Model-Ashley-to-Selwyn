# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:22:49 2017

@author: MichaelEK

Example script to extract vcsn data.
"""

from core.ecan_io.mssql import rd_squalarc
from pandas import read_csv, DataFrame, concat
from core.ecan_io.hilltop import rd_ht_wq_data

###########################################
### Parameters

#mtypes1 = ['Water Temperature', 'Ammonia Nitrogen']

sites_csv = r'E:\ecan\local\Projects\requests\michele\2017-10-02\sites.csv'
#poly = r'Y:\VirtualClimate\examples\example_area.shp'

mtype_names = ['temp', 'oxygen', 'turb', 'solids', 'nitrogen', 'nitrate', 'ammonium', 'phosp', 'coli', 'faecal', 'coliform', 'organic carbon']

from_date = '1997-07-01'
to_date = '2017-06-30'

convert_dtl = True
dtl_method = 'trend'

out_csv = r'E:\ecan\local\Projects\requests\michele\2017-10-02\wq_data_2017-10-04.csv'
out_sites_csv = r'E:\ecan\local\Projects\requests\michele\2017-10-02\wq_sites_2017-10-04.csv'

##########################################
### Get the data

sites = read_csv(sites_csv).site

## Get two parameters for the two sites during the time period and convert the dtl (taking half the dtl)
df2, sites_info = rd_ht_wq_data(sites=sites, convert_dtl=convert_dtl, dtl_method=dtl_method, output_site_data=True)


df3 = DataFrame()
sites_info2 = DataFrame()
for i in mtype_names:
    t3 = df2[df2['mtype'].str.contains(i, case=False)]
    s1 = sites_info[sites_info['mtype'].str.contains(i, case=False)]
    df3 = concat([df3, t3])
    sites_info2 = concat([sites_info2, s1])

t3 = df2[df2['mtype'].str.contains('pH')]
s1 = sites_info[sites_info['mtype'].str.contains('pH')]
df3 = concat([df3, t3])
sites_info2 = concat([sites_info2, s1])

#df2 = rd_squalarc(sites=sites, from_date=from_date, to_date=to_date, convert_dtl=convert_dtl, dtl_method=dtl_method)

#########################################
### Mangle

#df3 = df2.copy()
df3.loc[:, 'time'] = df3.loc[:, 'time'].dt.date
df4 = df3.pivot_table(values='data_dtl', index=['site', 'time'], columns='mtype', aggfunc='first')


#########################################
### Save the data

## Export the data as a csv file -- from the output dataframe
df4.to_csv(out_csv, encoding='utf-8', index=True)

sites_info2.to_csv(out_sites_csv, encoding='utf-8', index=False)





























