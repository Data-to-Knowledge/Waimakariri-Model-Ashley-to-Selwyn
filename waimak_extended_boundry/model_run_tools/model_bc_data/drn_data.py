# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 3:22 PM
"""
from __future__ import division
from waimak_extended_boundry.model_and_NSMC_build.m_packages.drn_packages import _get_drn_spd
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools import \
    mod_gns_model, get_model, get_race_data
import numpy as np
import os
import pickle
import pandas as pd
from env import sdp_required

def get_base_drn_spd(ncars=True):#todo test
    """
    get the drain data, but without the correct conductance values...
    :param ncars: boolean if True then return with teh north carpet drains(e.g. similar to the optimisation)
                  otherwise excluded these drains (used for some forawrd runs)
    :return:
    """
    data_path = os.path.join(sdp_required, 'base_drn_data.hdf')
    if ncars:
        drn_data = pd.read_hdf(data_path,'drn_with_n_car')
    else:
        drn_data = pd.read_hdf(data_path,'drn_without_n_car')
    return drn_data


def get_drn_spd(model_id): #todo move this up

    #todo returns record array with correct conductance values
    raise NotImplementedError

def get_drn_no_ncarpet_spd(model_id,recalc=False): # todo make this accessable

    #todo returns record array with correct conductance values
    pickle_path = "{}/model_{}_drn_wout_ncarpet.p".format(smt.temp_pickle_dir, model_id)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata
    model = get_model(model_id)
    org_drn_data = _get_drn_spd(1, 1)
    n_carpet = org_drn_data.loc[np.in1d(org_drn_data.group, ['ash_carpet', 'cust_carpet'])]
    n_carpet_ij = list(n_carpet.loc[:, ['i', 'j']].itertuples(False, None))
    outdata = model.get_package('drn').stress_period_data.data[0]
    idx = [e not in n_carpet_ij for e in zip(outdata['i'], outdata['j'])]

    outdata = outdata[idx]
    pickle.dump(outdata,open(pickle_path,'w'))

    return outdata

if __name__ == '__main__':
    print('done')
