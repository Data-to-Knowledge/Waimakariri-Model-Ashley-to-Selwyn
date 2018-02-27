# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 09:05:16 2017

@author: MichaelEK

Script to pull out data from Hydrotel.
"""

from core.ecan_io import rd_hydrotel

#########################################
#### Parameters

#sites = [66417, 66404]
sites = None
mtype = 'Rainfall'
#mtype = 'Flow'
from_date = '2017-07-21'
to_date = '2017-07-22 23:59'

resample = True
period = 'day'
n_periods = 2
fun = 'sum'
pivot = False

output_csv = r'S:\Surface Water\shared\projects\2017-07-23_flooding\precip\precip_data.csv'

########################################
#### Get data

data = rd_hydrotel(select=sites, mtype=mtype, from_date=from_date, to_date=to_date, resample=resample, period=period, n_periods=n_periods, fun=fun, pivot=pivot)


#######################################
#### Export data

data.to_csv(output_csv, header=True)




