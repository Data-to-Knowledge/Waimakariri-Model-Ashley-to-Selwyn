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


def get_base_drn_spd(ncars=True):
    """
    get the drain data, but without the correct conductance values which is model dependent
    metadata:
    the dataset of all drain BCs in the model domain including the northern carpet drains
        columns :
           'index': just an index no value
           'cond': the drain conductance, not a real value
           'elev': the drain invert elevation
           'group': the group for the drains
           'i': model row
           'j': model col
           'k': model layer
           'target_group': group for the target e.g. saltwater at factory rd
           'zone': which cwms zone
           'parameter_group': group for parameterisation some parameters without a target are grouped
                              with a reasonable set. e.g. for saltwater creek drain BCs below
                              last target's conductance is set to the nearest upstream target
           'mx': cell center x in the model
           'my': cell center y in the model

    :param ncars: boolean if True then return with teh north carpet drains(e.g. similar to the optimisation)
                  otherwise excluded these drains (used for some forawrd runs)
    :return:
    """
    data_path = os.path.join(sdp_required, 'base_drn_data.hdf')
    if ncars:
        drn_data = pd.read_hdf(data_path, 'drn_with_n_car')
    else:
        drn_data = pd.read_hdf(data_path, 'drn_without_n_car')
    return drn_data


def _get_drn_cond(model_id):
    """
    gets the drain conductance values for each drn parameter group for the given the model id
    :param model_id: identifier 'NsmcReal{nsmc_num:06d}'
    :return:
    """
    param_data = nc.Dataset(os.path.join(sdp_required, "nsmc_params_obs_metadata.nc"))

    nsmc_num = int(model_id[-6:])
    nidx = np.where(param_data.variables['nsmc_num'][:] == nsmc_num)[0][0]

    drns = np.array(param_data.variables['drns'][:])
    drn_cond = np.array(param_data.variables['drn_cond'][nidx, :])

    out = {k: v for k, v in zip(drns, drn_cond)}

    return out


def get_drn_spd(model_id, ncarpet=True, return_df=False):
    """
    get stress period data for the model, see get_base_drn for metadata
    :param model_id: identifier 'NsmcReal{nsmc_num:06d}'
    :param ncarpet: boolean return n carpet drains
    :param return_df: boolean if True return the pandas.DataFrame, otherwise return np.recarray,
                      which is set up to match the requirements of flopy stress period data
    :return:
    """

    drn_cond_mapper = _get_drn_cond(model_id)
    base_data = get_base_drn_spd(ncarpet)

    base_data.loc[:, 'parameter_group2'] = base_data.loc[:, 'parameter_group']
    base_data.replace({'parameter_group2': drn_cond_mapper}, inplace=True)
    base_data.loc[:, 'cond'] = base_data.loc[:, 'parameter_group2'].astype(float)
    base_data.drop('parameter_group2', axis=1, inplace=True)

    if return_df:
        return base_data
    else:
        return base_data.loc[:, ['k', 'i', 'j', 'elev', 'cond']].to_records(index=False)


def get_drn_no_ncarpet_spd(model_id):
    """
    depreciated here only to support backward comparability, correct access is through get_drn_spd
    :param model_id:
    :param recalc:
    :return:
    """
    return get_drn_spd(model_id, ncarpet=False)


if __name__ == '__main__':
    test = get_drn_spd('NsmcReal{:06d}'.format(17), False, True)

    print('done')
