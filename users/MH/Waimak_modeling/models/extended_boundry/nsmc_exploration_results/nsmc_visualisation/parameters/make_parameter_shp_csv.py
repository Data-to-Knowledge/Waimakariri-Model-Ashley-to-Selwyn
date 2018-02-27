# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 16/11/2017 8:51 AM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
import netCDF4 as nc
import os

def make_csvs_for_zeb(outdir):
    """
    a quiick script to give zeb data for picking rch and k pilot points to plot up together
    :param outdir:
    :return:
    """
    nc_file = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc"))

    rch_pt = np.array(nc_file.variables['rch_ppt'])
    rchx = np.array(nc_file.variables['rch_ppt_x'])
    rchy = np.array(nc_file.variables['rch_ppt_y'])
    rchgroup = np.array(nc_file.variables['rch_ppt_group'])
    meanings = nc_file.variables['rch_ppt_group'].flag_meanings.split(' ')
    values = nc_file.variables['rch_ppt_group'].flag_values
    rch_flag_mapper = {v:m for v,m in zip(values,meanings)}
    out_rch = pd.DataFrame({'id':rch_pt,'x':rchx,'y':rchy,'rch_group':rchgroup})
    out_rch = out_rch.replace({'rch_group':rch_flag_mapper})
    out_rch.to_csv(os.path.join(outdir,'rch_ppt.csv'))

    kv_pt = np.array(nc_file.variables['khv_ppt'])
    kv_x = np.array(nc_file.variables['khv_ppt_x'])
    kv_y = np.array(nc_file.variables['khv_ppt_y'])
    out_kv = pd.DataFrame({'id':kv_pt, 'x':kv_x, 'y':kv_y})
    out_kv.to_csv(os.path.join(outdir,'khv_ppt.csv'))

if __name__ == '__main__':
    make_csvs_for_zeb(r"T:\Temp\temp_gw_files")


