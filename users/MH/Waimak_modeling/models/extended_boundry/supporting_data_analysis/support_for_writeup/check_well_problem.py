# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/09/2018 9:25 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import numpy as np
from core.stats.LR_class import LR
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.interzone_n import get_chch_area_zones

if __name__ == '__main__':
    outpath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\model_checks\pumpin_problem"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    n_dataset = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\mednload_unc.nc")
    param_dataset = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")

    nsmc_nums = n_dataset.variables['nsmc_num']

    params = param_dataset.variables['pump_c'][np.in1d(param_dataset.variables['nsmc_num'][:],nsmc_nums)]

    chch_zones = get_chch_area_zones()

    all_n = n_dataset.variables['mednload'][:]

    ns = {}

    for layer in range(10):
        for zoneid, zone in chch_zones.items():
            temp = all_n[:,layer,zone] #todo check my indexing is rusty
            ns['zone {} layer {}'.format(zoneid,zone)] = temp.mean(axis=1).mean(axis=1)


    for id, val in ns.items():
        temp_lr = LR(params,val)
        temp_lr.plot(False,os.path.join(outpath,'{}.png'.format(id)))
    print('done')