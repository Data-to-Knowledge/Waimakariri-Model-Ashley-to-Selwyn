"""
Author: matth
Date Created: 14/06/2017 3:49 PM
"""

from __future__ import division
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import results_dir
import os
import numpy as np
import pandas as pd
import flopy

def output_k_conduct(version='m_strong_vert'):
    m = mt.get_base_mf_ss(model_version=version)

    base_dir = '{}/Material properties/ks_and_conductance_{}'.format(results_dir,version)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    kv = m.lpf.vka.array
    kh = m.lpf.hk.array
    for i in range(17):
        mt.array_to_raster('{}/model_kv_layer_{:02d}.tif'.format(base_dir,i+1),kv[i])
        mt.array_to_raster('{}/model_kh_layer_{:02d}.tif'.format(base_dir,i+1),kh[i])

    drain_conduct = np.zeros((190,365))
    str_conduct = np.zeros((190,365))
    str_data = pd.DataFrame(mt.get_base_str())
    for i in str_data.index:
        row, col = str_data.loc[i,['i','j']]
        str_conduct[row,col] += str_data.loc[i,'cond']

    str_conduct[np.isclose(str_conduct,0)] = -99
    drns = pd.DataFrame(m.drn.stress_period_data.data[0])

    for i in drns.index:
        row, col = drns.loc[i, ['i', 'j']]
        drain_conduct [row, col] += drns.loc[i,'cond']

    drain_conduct[~mt.get_drn_samp_pts_dict()['all_drains']] = -99

    mt.array_to_raster('{}/drain_conduct.tif'.format(base_dir), drain_conduct)
    mt.array_to_raster('{}/str_conduct.tif'.format(base_dir), str_conduct)

def output_heads_mfss(version='m_strong_vert'):
    base_dir = "{}/Material properties/heads_ss_{}".format(results_dir,version)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    hds = flopy.utils.HeadFile(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\base_model_runs\base_ss_mf\base_SS.hds")
    head_array = hds.get_data(kstpkper=(0,0))
    for i in range(17):
        mt.array_to_raster('{}/model_heads_SS_layer_{:02d}.tif'.format(base_dir, i + 1), head_array[i])

if __name__ == '__main__':
    #output_k_conduct()
    output_heads_mfss()

