# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:22:49 2017

@author: MichaelEK

Example script to extract vcsn data.
"""

from core.ecan_io.met import rd_niwa_vcsn
from numpy import nan, full, arange, seterr
seterr(invalid='ignore')

###########################################
### Parameters

mtypes2 = ['precip', 'PET']

sites = ['P146097']
poly = r"P:\Hurunui Waiau SRP\St Anne's Lagoon"

include_sites = True

initial_storage = 800
max_storage = 830
min_storage = 0

out_csv = r"P:\Hurunui Waiau SRP\St Anne's Lagoon\bucket_model_results1.csv"

##########################################
### Get the data

## Get the precip data for the two sites, include the sites in the output, and convert the x and y coordinates from wgs 84 to NZTM
df1 = rd_niwa_vcsn(mtypes=mtypes2, sites=sites, include_sites=include_sites)

##########################################
### Bucket model

time = df1.time
rain = df1.rain.values
et = df1.pe.values

storage = full(len(time), nan)
storage[-1] = initial_storage

for i in arange(len(time)):
    storage[i] = rain[i] - et[i] + storage[i - 1]
    if storage[i] > max_storage:
        storage[i] = max_storage
    if storage[i] < min_storage:
        storage[i] = min_storage

df2 = df1[['time', 'rain', 'pe']].copy()
df2.loc[:, 'storage'] = storage
df2.set_index('time', inplace=True)

#########################################
### Plot

ax = df2[['rain', 'pe']].plot(alpha=0.5)
df2[['storage']].plot(ax=ax, secondary_y=True, style='black')



#########################################
### Save the data

## Export the data as a csv file
df2.to_csv(out_csv, index=False)






























