# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 16:07:49 2017

@author: MichaelEK
"""

from pandas import read_csv, to_datetime, DateOffset, concat


#################################################
#### Input parameters

lsrm_csv = r'E:\ecan\shared\projects\lsrm\output_waipara1.csv'
ds_csv = r'E:\ecan\shared\projects\lsrm\waipara_ds_test.csv'

###############################################
#### Read in data

lsrm_results = read_csv(lsrm_csv)
ds_results = read_csv(ds_csv)

ds_results.loc[:, 'date'] = to_datetime(ds_results.date, format='%Y%m') + DateOffset(months=1) - DateOffset(days=1)
ds_results.set_index('date', inplace=True)

lsrm_results.rename(columns={'time': 'date'}, inplace=True)
lsrm_results.loc[:, 'date'] = to_datetime(lsrm_results.loc[:, 'date'])
lsrm_results.set_index('date', inplace=True)

############################################
#### Aggregate lsrm results

lsrm1 = lsrm_results[['precip', 'pet', 'site_area']].copy()
lsrm1.loc[:, 'aet2'] = lsrm_results['irr_aet'] + lsrm_results['non_irr_aet']
lsrm1.loc[:, 'drainage2'] = lsrm_results['irr_drainage'] + lsrm_results['non_irr_drainage']
lsrm1.loc[:, 'demand2'] = lsrm_results['irr_demand']

lsrm2 = lsrm1.groupby(level='date').sum()
lsrm3 = lsrm2.div(lsrm2.site_area, axis=0).drop('site_area', axis=1)
lsrm4 = lsrm3.resample('M').sum()

lsrm5 = lsrm4.rename(columns={'precip': 'precip2', 'pet': 'pet2', 'irr_demand': 'demand2'}).round(2)
lsrm5[lsrm5.isnull()] = 0

############################################
#### Combine data sets

both1 = concat([ds_results, lsrm5], axis=1, join='inner')
both2 = concat([ds_results, lsrm5], axis=1, join='inner')




###########################################
#### Plot

both1[['demand', 'demand2']].plot()
both1[['drain', 'drainage2']].plot()

both2[['demand', 'demand2']].plot()
both2[['drain', 'drainage2']].plot()




#######################################
#### Testing

p1 = new_rain_et.copy()
p1.index = p1.index.droplevel(['x', 'y'])
p2 = p1.groupby(level='time').mean()
p3 = p2.resample('M').sum()

m1 = model_var.copy()
m1.set_index('time', inplace=True)
m2 = m1.resample('M')[['precip', 'pet', 'site_area']].sum()
m2.div(m2.site_area, axis=0)

p1 = input1[['time', 'precip', 'pet', 'site_area']].copy()
p1.set_index('time', inplace=True)
p3 = p1.resample('M').sum()
p3.div(p3.site_area, axis=0)

p1 = new_rain_et1.copy()
p1.index = p1.index.droplevel(['x', 'y'])
p2 = p1.groupby(level='time').mean()
p3 = p2.resample('M').sum()


p1 = lsrm_results.copy()
p3 = p1.resample('D').sum()
p3.div(p3.site_area, axis=0)
































