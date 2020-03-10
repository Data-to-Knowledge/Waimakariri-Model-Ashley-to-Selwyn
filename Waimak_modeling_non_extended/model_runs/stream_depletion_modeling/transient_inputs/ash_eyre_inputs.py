"""
Author: matth
Date Created: 22/05/2017 10:25 AM
"""

from __future__ import division
from core.classes.hydro import hydro
import pandas as pd
import numpy as np
import os
import pickle
from users.MH.Waimak_modeling.supporting_data_path import sdp

def get_monthly_str_scale(recalc=False):
    picklepath = '{}/inputs/pickeled_files/mean_year_stream_scalings.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        str_scaling = pickle.load(open(picklepath))
        return str_scaling


    site_to_num_dict = {
        'ashley': 66204,
        'eyre': 166405,
        'a_tribs': 66213
    }

    h1 = hydro().get_data(mtypes='flow', sites=site_to_num_dict.values(), from_date='01/01/2008', to_date='31/12/2016')
    str_scaling = {}
    for site in site_to_num_dict.keys():
        id_num = site_to_num_dict[site]
        temp = pd.DataFrame(h1.data.loc['flow',id_num])
        temp['year'] = temp.index.year
        temp['month'] = temp.index.month
        data_mean = temp.data.mean()
        g = pd.groupby(temp,['year','month'])
        temp2 = g.aggregate({'data': np.mean}).reset_index()
        temp2.data *= 1/data_mean
        g2 = pd.groupby(temp2,['month'])
        temp3 = g2.aggregate({'data': np.mean})
        temp3.data *= 12/temp3.data.sum() # scale so no flow lost
        str_scaling[site] = temp3.data

    pickle.dump(str_scaling,open(picklepath,'w'))
    return str_scaling

if __name__ == '__main__':
    test = get_monthly_str_scale()
    print(test)
