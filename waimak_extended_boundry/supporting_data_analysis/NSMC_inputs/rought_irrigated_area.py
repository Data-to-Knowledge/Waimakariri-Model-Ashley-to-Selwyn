# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/08/2017 10:17 AM
"""

from __future__ import division
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_and_NSMC_build.m_packages.rch_packages import _get_rch
import numpy as np


no_flow = smt.get_no_flow(0)
no_flow[no_flow<0]=0
zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
zones[~no_flow.astype(bool)] = np.nan
# waimak = 4, chch_wm = 7, selwyn=8
w_idx = np.isclose(zones, 4)
c_idx = np.isclose(zones, 7)
s_idx = np.isclose(zones, 8)
irrigation = smt.shape_file_to_model_array(r"C:\Users\MattH\Downloads\irrigated_Area.shp",'OBJECTID',True)
irrigation[~no_flow.astype(bool)] = np.nan
confined = smt.shape_file_to_model_array(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\confined_zone_full_extent.shp",'Id',True)
confined[~no_flow.astype(bool)] = np.nan
rch = _get_rch()
rch[~no_flow.astype(bool)] = np.nan
for idx, name in zip([w_idx,c_idx,s_idx],['waimak','chch','selwyn']):
    print('')
    print(name)
    irrigated_area = ((np.isfinite(irrigation[idx])) & (~np.isfinite(confined[idx]))).sum()
    total_area = idx.sum()
    conf = np.isfinite(confined[idx]).sum()
    non_irr = total_area - irrigated_area - conf
    print('total_area {} cells'.format(total_area))
    print('irrgated {}% {} cells'.format(irrigated_area/total_area,irrigated_area))
    print('non-irrgated {}% {} cells'.format(non_irr/total_area,non_irr))
    print('confined {}% {} cells'.format(conf/total_area,conf))
    print((conf+non_irr+irrigated_area)/total_area)
    print('average recharge: {}'.format(rch[idx].mean()*1000*365))

dryland =rch[np.isnan(confined) & np.isnan(irrigation) & (no_flow.astype(bool))].mean()*1000*365
print ('average rch in dryland: {} mm/yr'.format(dryland))
confined_drl = rch[np.isfinite(confined) & np.isnan(irrigation)].mean()*1000*365
print ('average rch in confined dryland: {} mm/yr'.format(confined_drl))
waiir = rch[np.isfinite(irrigation) & (zones==4)].mean()*1000*365
print ('average rch in waimak irrigated: {} mm/yr'.format(waiir))
sel_ir = rch[np.isfinite(irrigation) & (~(zones==4))].mean()*1000*365
print ('average rch in selwyn irrigated: {} mm/yr'.format(sel_ir))
