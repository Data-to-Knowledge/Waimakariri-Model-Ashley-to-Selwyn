# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/08/2017 10:17 AM
"""

from __future__ import division
from waimak_extended_boundry import smt
import numpy as np

def get_rch_index_array():
    no_flow = smt.get_no_flow(0)
    no_flow[no_flow<0]=0
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    zones[~no_flow.astype(bool)] = np.nan
    # waimak = 4, chch_wm = 7, selwyn=8
    irrigation = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\irrigated_Area.shp".format(smt.sdp),'OBJECTID',True)
    irrigation[~no_flow.astype(bool)] = np.nan
    confined = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\confined_zone_full_extent.shp".format(smt.sdp),'Id',True)
    confined[~no_flow.astype(bool)] = np.nan

    # 0=ibound 1 = dryland, 2=confined dryland, 3= selwyn irrigated, 4 =waimak irrigated
    out_array = np.ones((smt.rows,smt.cols))
    out_array[~no_flow.astype(bool)] = 0
    out_array[np.isfinite(confined)] =  2
    out_array[np.isfinite(irrigation) & ~np.isclose(zones,4)] =  3
    out_array[np.isfinite(irrigation) & np.isclose(zones,4)] =  4
    return out_array



