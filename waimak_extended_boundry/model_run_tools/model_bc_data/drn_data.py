# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 3:22 PM
"""
from __future__ import division
import numpy as np
import os
import pandas as pd
from env import sdp_required
import netCDF4 as nc

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

def _get_drn_cond(model_id):
    """
    gets the drain conductance values for the model id
    :param model_id: identifier 'NsmcReal{nsmc_num:06d}'
    :return:
    """
    param_data = nc.Dataset(os.path.join(sdp_required,"nsmc_params_obs_metadata.nc"))

    nsmc_num = int(model_id[-6:])
    nidx = np.where(param_data.variables['nsmc_num'][:] == nsmc_num)[0][0]

    drns = np.array(param_data.variables['drns'][:])
    drn_cond = np.array(param_data.variables['drn_cond'][nidx,:])

    out = {k:v for k,v in zip(drns,drn_cond)}

    return out

def get_drn_spd(model_id, ncarpet=True): #todo test
    """
    get stress period data for the model with n carpert drains
    :param model_id: identifier 'NsmcReal{nsmc_num:06d}'
    :param ncarpet: boolean return n carpet drains
    :return:
    """

    drn_cond_mapper = _get_drn_cond(model_id)
    base_data = get_base_drn_spd(ncarpet)

    base_data.replace({'parameter_group':drn_cond_mapper})
    base_data.loc[:,'cond'] = base_data.loc[:, 'parameter_group'].astype(float)

    return base_data.loc[:,['k','i','j','elev','cond']].to_records()


def get_drn_no_ncarpet_spd(model_id):
    """
    this is here only to support backward compatablity
    :param model_id:
    :param recalc:
    :return:
    """
    return get_drn_spd(model_id,ncarpet=False)

if __name__ == '__main__':
    get_drn_spd('NsmcReal{:06d}'.format(15))

    print('done')
