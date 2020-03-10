# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/08/2017 4:11 PM
"""

from __future__ import division
import numpy as np
from waimak_extended_boundry import smt
from waimak_extended_boundry import _get_rch
import matplotlib.pyplot as plt

new_no_flow = smt.get_no_flow()
zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
zones[~new_no_flow[0].astype(bool)] = 0
# waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
confin_rch_zone = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/chch_wm_rch_split_chch_form.shp".format(smt.sdp),
                                                'ID', True)
zones[(zones == 7) & (np.isfinite(confin_rch_zone))] = 9

model_rch = _get_rch()*365*1000
ds_rch = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp),'rch',True)*365
ds_paw = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp),'paw',True)
print('full infultration with no PAW')
for zone, name in zip([4,7,9,8],['waimak','chch_wm','chch_form','selwyn']):
    idx = np.isclose(zones.astype(float),zone)
    print('{} zone: model:{}, david scott:{}'.format(name,model_rch[idx].mean(),ds_rch[idx].mean()))

ds_rch[ds_paw<=400000] = 0
print('no influtration with PAW')
for zone, name in zip([4,7,9,8],['waimak','chch_wm','chch_form','selwyn']):
    idx = np.isclose(zones.astype(float),zone)
    print('{} zone: model:{}, david scott:{}'.format(name,model_rch[idx].mean(),ds_rch[idx].mean()))
fig, (ax,ax2) = plt.subplots(2)
smt.plt_matrix(model_rch,ax=ax)
smt.plt_matrix(ds_rch,ax=ax2)
plt.show(fig)

