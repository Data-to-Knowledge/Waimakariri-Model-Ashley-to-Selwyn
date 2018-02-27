"""
Author: matth
Date Created: 17/03/2017 2:26 PM
"""

from __future__ import division
import pandas as pd
import numpy as np
from users.MH.Well_classes import qual_parameters


data = pd.read_csv(r"P:\Groundwater\Annual groundwater quality survey 2016\data\all_output.csv")
data = data[data.NO3_N_trend.notnull()]

# move the below into data processing and not getting exactly the same results as hill's

# charge check
charge_balance_pos = np.zeros((len(data.index)))
charge_balance_neg = np.zeros((len(data.index)))
for param in qual_parameters:
    if qual_parameters[param]['mol_mass'] is None:
        continue
    temp = data[param]/qual_parameters[param]['mol_mass']*qual_parameters[param]['charge']
    if qual_parameters[param]['charge'] < 0:
        charge_balance_neg += temp
    else:
        charge_balance_pos += temp

data['charge_bal'] = (charge_balance_pos + charge_balance_neg)/((charge_balance_pos - charge_balance_neg))


numwells = len(data.depth)
depth = data.depth
print('{} wells less than 20 '.format(len(depth[depth<=20])))
print('{} percentage wells less than 20m'.format(len(depth[depth<=20])/len(data.depth)))
print('{} wells 20-50m deep '.format(len(depth[(depth>20)&(depth<=50)])))
print('{} wells less than 500-100m deep'.format(len(depth[(depth>50)&(depth<=100)])))
print('{} wells >100'.format(len(depth[(depth>50)&(depth<=100)])))
print('maximum depth {}'.format(depth.max())) #may need to do more reserach on this one if the deep well is ever dropped from sampling

numtrendwells=len(data.NO3_N_trend[data.NO3_N_trend.isin(['increasing','decreasing','no trend'])])
inc = len(data.NO3_N_trend[data.NO3_N_trend == 'increasing'])
dec = len(data.NO3_N_trend[data.NO3_N_trend == 'decreasing'])
notrend = len(data.NO3_N_trend[data.NO3_N_trend == 'no trend'])
print('increasing N {} {}%'.format(inc,inc/numtrendwells))
print('decreasing  N {} {}%'.format(dec,dec/numtrendwells))
print('no trend N {} {}%'.format(notrend,notrend/numtrendwells))

for val in set(data.cwms_zone):
    print('NO3 trend by zone')
    print(val)
    print (len(data.NO3_N_trend[(data.NO3_N_trend == 'increasing') & (data.cwms_zone == val)]))

numtrendwells=len(data.drp_trend[data.drp_trend.isin(['increasing','decreasing','no trend'])])
inc = len(data.drp_trend[data.drp_trend == 'increasing'])
dec = len(data.drp_trend[data.drp_trend == 'decreasing'])
notrend = len(data.drp_trend[data.drp_trend == 'no trend'])
print('increasing P {} {}%'.format(inc,inc/numtrendwells))
print('decreasing  P {} {}%'.format(dec,dec/numtrendwells))
print('no trend P {} {}%'.format(notrend,notrend/numtrendwells))

for val in set(data.cwms_zone):
    print('P trend by zone')
    print(val)
    print (len(data.drp_trend[(data.drp_trend == 'increasing') & (data.cwms_zone == val)]))

print('num wells {} percentage {} with ecoli detection'.format(len(data.ecoli[data.ecoli > 0]),len(data.ecoli[data.ecoli > 0])/numwells))
print('num wells {} percentage {} with ecoli detection in wells '
      'less than 20m'.format(len(data.ecoli[(data.ecoli > 0) & (data.depth < 20)]),len(data.ecoli[(data.ecoli > 0) & (data.depth < 20)])/numwells))
print('num wells {} ecoli detection wells over 20'.format(len(data.ecoli[(data.ecoli > 0) & (data.depth >=20)])))


print ('done')