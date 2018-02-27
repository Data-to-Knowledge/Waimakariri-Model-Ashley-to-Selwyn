# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/08/2017 12:07 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.rch_packages import _get_rch
from copy import deepcopy


if __name__ == '__main__':
    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp),'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = np.nan
    confin_rch_zone = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/chch_wm_rch_split_chch_form.shp".format(smt.sdp),'ID',True)
    zones[(zones==7) & (np.isfinite(confin_rch_zone))] = 9
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    w_idx = np.isclose(zones,4)
    c_idx = np.isclose(zones,7)
    s_idx = np.isclose(zones,8)
    cc_idx = np.isclose(zones,9)
    all_idx = np.isfinite(zones)

    outdata = pd.DataFrame(columns=['waimak','selwyn','chch_wm','chch_wm_con', 'total'])
    rch = _get_rch()
    irrg = np.isfinite(smt.shape_file_to_model_array(r"T:\Temp\temp_gw_files\Irrigated_land_2016.shp",'OBJECTID',True))

    irrigation = deepcopy(rch)
    irrigation[~irrg] = np.nan
    dryland = deepcopy(rch)
    dryland[irrg] = np.nan

    rch_ds = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp), 'rch', True)/1000


    irrigation_ds = deepcopy(rch_ds)
    irrigation_ds[~irrg] = np.nan
    dryland_ds = deepcopy(rch_ds)
    dryland_ds[irrg] = np.nan
    for name,dat in zip(['irrigation','dryland', 'irrigation_ds','dryland_ds'], [irrigation, dryland,irrigation_ds,dryland_ds]):
        for idx, zone in zip([w_idx,c_idx,s_idx,cc_idx,all_idx],['waimak','chch_wm', 'selwyn', 'chch_wm_con', 'total']):
            outdata.loc[name,zone] = np.nanmean(dat[idx])

    print outdata*365*1000 #convert to mm/year