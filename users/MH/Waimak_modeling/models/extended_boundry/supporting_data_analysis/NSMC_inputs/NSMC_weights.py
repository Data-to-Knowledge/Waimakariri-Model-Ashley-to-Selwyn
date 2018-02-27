# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 1/09/2017 3:17 PM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.targets.gen_target_arrays import get_vertical_gradient_targets, get_head_targets
from users.MH.Waimak_modeling.models.extended_boundry.targets.briochs_vert_tar_script import get_vert_targets_full
import numpy as np

target_sd_actual = { # in m3/s
    #sfr flow
'sfo_c_benn':0.06,
'sfo_c_oxf':0.1458,
'sfo_c_swan':0.1275,
'sfo_c_tlks':0.0875,
'sfo_c_skew': 0.0875,
'sfo_e_wolf': 0.1,
'sfo_e_poyn': 0.05,
'sfo_e_down': 0.005,
'sfo_e_seyr': 0.0005,

#SFR flux
'sfx_e1_stf':0, # not including as it is dependent on SFR influx
'sfx_w1_cou':1.36,
'sfx_w2_tom':1.36,
'sfx_w3_ros':1.36,
'sfx_w4_mcl':1.36,
'sfx_w5_wat':1.36,
'sfx_w6_sh1':1.36,
'sfx_w_all':3.33,
'sfx_e_all':0, # not including as it is dependent on SFR influx which we are varing

#drn flux
'd_ash_est': 0.075,
'd_salt_top': 0.0675,
'd_salt_fct':0.035,
'd_kuku_leg':0.1125,
'd_tar_gre':0.0275,
'd_tar_stok':0.025,
'd_cam_mrsh':0.0425,
'd_nbk_mrsh':0.165,
'd_sbk_mrsh':0.0425,
'd_cam_yng':0.102247,
'd_cam_revl':0.1,
'd_oho_mlbk':0.0375,
'd_oho_whit':0.0025,
'd_oho_jefs':0.0075,
'd_oho_misc':0.05,
'd_oho_btch':0.0725,
'd_oho_kpoi':0.0475,
'd_emd_gard':0.045,
'd_sil_harp':0.018,
'd_sil_heyw':0.0975,
'd_sil_ilnd':0.2225,
'd_cour_nrd':0.0925,
'd_smiths':0.03,
'd_kairaki':0.04,


#chb flux
'chb_ash':0, # not including as too uncertain
'chb_cust':0.162,
'chb_chch':0.383,
'sel_off': 0, # not including as too uncertain

#calculated from relative
#sfr flux
'mid_ash_g': 5.20 * 0.29, # 29%
'sfx_a4_sh1': 0.4 * 0.29, # 29%

#drain flux
'sel_str': 9.8 * 0.10, # 10%
'chch_str': 10 * 0.125 # 12.5%
}

def get_NSMC_weights():
    actual = pd.DataFrame(data=target_sd_actual.values(),index=target_sd_actual.keys(),columns=['std'])
    actual.loc[:,'std'] *=86400 # convert from m3/s to m3/day
    actual.loc[:,'weight'] = 1/actual.loc[:,'std']
    actual.loc[~np.isfinite(actual.weight),'weight'] = 0



    #vertical gradient targets
    vert_targets = get_vert_targets_full()
    vert_targets = pd.DataFrame(vert_targets.set_index('name').loc[:,'cwgt'])
    vert_targets = vert_targets.rename(columns={'cwgt': 'weight'})
    # head targets
    heads = get_head_targets()
    heads = pd.DataFrame(heads.weight)

    outdata = pd.concat((actual,vert_targets,heads))
    outdata = outdata.loc[:,'weight']
    return outdata

if __name__ == '__main__':
    test = get_NSMC_weights()
    print test