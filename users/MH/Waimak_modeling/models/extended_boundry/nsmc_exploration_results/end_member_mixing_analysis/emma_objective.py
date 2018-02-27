# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/11/2017 2:20 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd


def calc_objective_function(modeled, target, weights, groups, group_weights):
    """

    :param modeled: array of modeled values (nsmc, value)
    :param target: array of target values (value,)
    :param weights: array of un-altered weights (value)
    :param groups: array of group ids or None (no group weighting) (value)
    :param group_weights: dictionary of group ids: multiplicitive weight or None no group weighting
    :return:
    """
    # some checks

    if not (groups is None and group_weights is None) and not (groups is not None and group_weights is not None):
        raise ValueError('groups and groupweights must both be None or Not none')

    if set(groups) != set(group_weights.keys()):
        raise ValueError('all group ids must be in keys')

    residuals = modeled - target[np.newaxis, :]

    if groups is not None:
        weights = np.array([w * group_weights[g] for w, g in zip(weights, groups)])
    weights = weights[np.newaxis, :]
    phi = np.power(weights * residuals, 2,).sum(axis=1)

    return phi


def create_emma_phis(outpath):
    data = pd.read_csv(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\nsmc_results\emma_filter\river_end_member_mixing.csv"),
        index_col=0)

    target_data = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\4EM_WaimakResults_NSMCFilterWithStreams.csv",
        index_col=0)
    target_data.loc[:, 'Group'] = [e.lower() for e in target_data.Group]
    target_data.loc[:, 'SuperGroup'] = [e.lower() for e in target_data.SuperGroup]
    # river_5th_4EM	river_med_4EM	river_95th_4EM

    target_data.loc[:, 'sd'] = ((target_data.loc[:, 'river_med_4EM'] - target_data.loc[:, 'river_5th_4EM']) +
                                (target_data.loc[:, 'river_95th_4EM'] - target_data.loc[:, 'river_med_4EM'])) / 2

    targets = target_data.groupby('Group').aggregate({'sd': np.mean, 'river_med_4EM': np.mean})
    targets.loc[:, 'supergroup'] = [target_data.loc[target_data.Group == g, 'SuperGroup'].iloc[0] for g in
                                    targets.index]
    sites = list(set(target_data.Group).intersection(data.keys()) - {'inland'})
    weights = 1 / targets.loc[sites, 'sd'].values

    #{'deep eyreton kaiapoi belfast',
    # 'deep mcleans avonhead',
    # 'eyrewell forest',
    # 'mid mid chch',
    # 'mid wm airfield',
    # 'shallow mid chch',
    # 'shallow n chch',
    # 'streams'}
    # weight up by 2 orders of magnitutde
    group_weights = {'no_weighting': _make_weight_dict([]),
                     'equal_num': _make_weight_dict('equal',targets),
                     'chch_weighted': _make_weight_dict(['shallow n chch', 'shallow mid chch', 'mid mid chch',
                                                         'deep mcleans avonhead']),
                     'stream_weighted': _make_weight_dict(['streams']),
                     'ewf_weighted': _make_weight_dict(['eyrewell forest']),
                     }

    outdata = {}
    for key, g_wght in group_weights.items():
        temp_phi = calc_objective_function(modeled=data.loc[:, sites].values, target=targets.loc[sites,'river_med_4EM'].values, weights=weights,
                                groups=targets.loc[sites, 'supergroup'].values, group_weights=g_wght)
        outdata[key] = temp_phi
    outdata = pd.DataFrame(data=outdata,index=data.index)
    outdata.to_csv(outpath)


def _make_weight_dict(keys, targets=None):
    base_dict = {'deep eyreton kaiapoi belfast': 1,
                 'deep mcleans avonhead': 1,
                 'eyrewell forest': 1,
                 'mid mid chch': 1,
                 'mid wm airfield': 1,
                 'shallow mid chch': 1,
                 'shallow n chch': 1,
                 'streams': 1}

    if keys == 'equal':
        nums = targets.groupby('supergroup').count()
        for key in base_dict:
            base_dict[key] = 1/nums.loc[key,'sd'] #sd is just a filler from the groupby, I'm being lazy... or panicked

    else:
        for key in keys:
            base_dict[key] = 100  # raise weighting by two orders of magnitude

    return base_dict

def make_n_converged():
    import socket
    import netCDF4 as nc
    if socket.gethostname() != 'GWATER02':
        raise ValueError('must be run on GWater02')

    nc_data = nc.Dataset(r"C:\mh_waimak_model_data\mednload_ucn.nc")
    conv = []
    for i, num in enumerate(np.array(nc_data.variables['nsmc_num'])):
        if not np.isnan(np.array(nc_data.variables['mednload'][i,1])).all():
            conv.append(num)

    with open(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\nsmc\mednload_converged.txt"),'w') as f:
        f.writelines(['{}\n'.format(e) for e in conv])



if __name__ == '__main__':
    make_n_converged()