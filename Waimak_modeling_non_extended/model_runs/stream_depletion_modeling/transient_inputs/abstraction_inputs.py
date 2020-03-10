"""
Author: matth
Date Created: 22/05/2017 10:26 AM
"""

from __future__ import division
from core.classes.hydro import hydro
import users.MH.Waimak_modeling.model_tools as mt
import pandas as pd
import numpy as np
import pickle
import os
from users.MH.Waimak_modeling.supporting_data_path import sdp


def get_mean_month_well_scaling(recalc=False):
    picklepath = '{}/inputs/pickeled_files/mean_year_abstraction_scalings.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        use_scaling = pickle.load(open(picklepath))
        return use_scaling

    wells = mt.get_all_well_data()
    sites = list(wells.index)


    h1 = hydro().get_data(mtypes='usage', sites=sites)

    outdata = pd.DataFrame(index=wells.index)
    outdata['type'] = wells['type']
    outdata['use1'] = wells['use1']
    outdata['use2'] = wells['use2']
    outdata[1] = np.nan
    outdata[2] = np.nan
    outdata[3] = np.nan
    outdata[4] = np.nan
    outdata[5] = np.nan
    outdata[6] = np.nan
    outdata[7] = np.nan
    outdata[8] = np.nan
    outdata[9] = np.nan
    outdata[10] = np.nan
    outdata[11] = np.nan
    outdata[12] = np.nan

    null_count = 0
    for well in wells.index:
        if well not in h1.data.index.get_level_values('site'):
            continue

        temp = pd.DataFrame(h1.data.loc['usage',well])
        data_mean = temp.data.mean()
        if data_mean == 0:
            continue
        temp['month'] = temp.index.month
        temp['year'] = temp.index.year
        g = pd.groupby(temp, ['year', 'month'])
        temp2 = g.aggregate({'data': np.mean}).reset_index()
        temp2.data *= 1/data_mean
        g2 = pd.groupby(temp2,['month'])
        temp3 = g2.aggregate({'data':np.mean})
        if len(temp3) != 12:
            null_count += 1
            continue

        outdata.loc[well,[1,2,3,4,5,6,7,8,9,10,11,12]] = list(temp3.data)

    uses = ['irrigation', 'public_water_supply']
    use_scaling = {}

    for use in uses:
        temp = outdata[(outdata['use1'] == use) | (outdata['use2'] == use)]
        use_scaling[use] = temp.mean(axis=0)*12/temp.mean(axis=0).sum()

    pickle.dump(use_scaling,open(picklepath,'w'))

    return use_scaling

if __name__ == '__main__':
    get_mean_month_well_scaling(True)