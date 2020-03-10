"""
Author: matth
Date Created: 28/04/2017 9:45 AM
"""

from __future__ import division
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.supporting_data_path import sdp
from core.ecan_io import rd_sql, sql_db
# test if we can simply re-build the wells information

well_details = rd_sql(**sql_db.wells_db.well_details)

all_consent_data = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp)) # allocation data from mike

well_s_waimak_data = pd.read_csv("{}/inputs/wells/consents_s_waimak.csv".format(sdp)) # consent data from Brioch

# figure out how many wells per consent

consents = {}
for con in well_s_waimak_data.consent:
    wells = all_consent_data.wap[(all_consent_data.crc == con) & (all_consent_data.take_type=='Take Groundwater')]
    consents[con] = wells

num_wells = np.array([len(e) for e in consents.values()])

print('max: {}'.format(num_wells.max()))
print('mean: {}'.format(num_wells.mean()))
print('min: {}'.format(num_wells.min()))
for i in range(0,num_wells.max()+1):
    print('{} wells: {}'.format(i, len(num_wells[num_wells==i])))


for key in consents.keys():
    if len(consents[key]) in [0,8]:
        print("{}: {}".format(key,consents[key]))

distance = {}
for i in well_s_waimak_data.index:
    con = well_s_waimak_data.loc[i,'consent']
    if len(consents[con]) == 0:
        continue

    well = np.atleast_1d(consents[con])[0] #for now just look at the first well
    con_x = well_s_waimak_data.loc[i,'x']
    con_y = well_s_waimak_data.loc[i,'y']
    well_x = well_details.NZTMX[well_details.WELL_NO == well].iloc[0]
    well_y = well_details.NZTMY[well_details.WELL_NO == well].iloc[0]

    distance[con] = ((well_x-con_x)**2 + (well_y - con_y)**2)**0.5

dists = np.array(distance.values())
print('max_dis: {}'.format(dists.max()))
print('mean_dis: {}'.format(dists.mean()))
print('min_dis: {}'.format(dists.min()))
print('num over 1km: {}'.format(len(dists[dists > 1000])))
print(dists[dists > 1000])

print('done')


