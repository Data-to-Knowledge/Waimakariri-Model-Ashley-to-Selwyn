"""
Author: matth
Date Created: 29/05/2017 9:09 AM
"""

from __future__ import division
from core.classes.hydro import hydro
import pandas as pd
import numpy as np
import pickle
import os
from users.MH.Waimak_modeling.supporting_data_path import sdp

def get_s_cal_data(recalc=False):

    picklepath = '{}/inputs/pickeled_files/s_calibration_targets.p'.format(sdp)

    if os.path.exists(picklepath) and not recalc:
        s_cal_targets = pickle.load(open(picklepath))
        return s_cal_targets


    site_to_num_dict = {
        'kaiapoi_harpers': 66415,
        'cust_threlkelds': 66417,
        'cam_youngs': 66409
    }

    # get data
    h1 = hydro().get_data(mtypes='flow', sites=site_to_num_dict.values(), from_date='01/01/2008', to_date='31/12/2016')
    s_cal_targets = pd.DataFrame()
    for site in site_to_num_dict.keys():
        temp = pd.DataFrame(h1.data.loc['flow',site_to_num_dict[site],:])
        temp['month'] = temp.index.month
        g = pd.groupby(temp,['month'])
        temp2 = g.aggregate({'data': np.mean})
        s_cal_targets[site] = temp2.loc[:,'data'] * 86400


    pickle.dump(s_cal_targets,open(picklepath,'w'))
    return s_cal_targets



if __name__ == '__main__':
    s_cal = get_s_cal_data(recalc=True)