# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/09/2018 9:25 AM
"""

from __future__ import division
import env
import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from scipy.stats import spearmanr, pearsonr
from waimak_extended_boundry.model_runs.n_analysis.interzone_n import get_chch_area_zones

if __name__ == '__main__':
    outpath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\model_checks\pumpin_problem"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    n_dataset = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\mednload_unc.nc"))
    param_dataset = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc"))

    nsmc_nums = n_dataset.variables['nsmc_num']

    params = param_dataset.variables['pump_c'][np.in1d(param_dataset.variables['nsmc_num'][:],nsmc_nums)]

    chch_zones = get_chch_area_zones()

    all_n = np.array(n_dataset.variables['mednload'])

    ns = {}
    outdata = pd.DataFrame(columns=['pearsonr','pearsonp','spearmanr', 'spearmanp'])
    for layer in range(10):
        for zoneid, zone in chch_zones.items():
            temp = all_n[:,layer,zone]
            ns['layer_{:02d}_zone_{}'.format(layer, zoneid)] = np.nanmean(temp, axis=1)

    for id, val in ns.items():
        idx = np.isfinite(val)
        outdata.loc[id,'pearsonr'], outdata.loc[id,'pearsonp'] = pearsonr(params[idx],val[idx])
        temp = spearmanr(params[idx],val[idx])
        outdata.loc[id,'spearmanr'] = temp.correlation
        outdata.loc[id,'spearmanp'] = temp.pvalue
        if False:
            temp_lr = LR(params[idx],val[idx])
            temp_lr.plot(False,os.path.join(outpath,'{}.png'.format(id)))
    outdata.to_csv(os.path.join(outpath,'correlations.csv'))
    print('done')