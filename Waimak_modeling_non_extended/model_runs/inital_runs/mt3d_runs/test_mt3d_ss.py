"""
Author: matth
Date Created: 16/05/2017 8:22 AM
"""

from __future__ import division
import numpy as np
from copy import deepcopy
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, sdp


base_dir = '{}/component_con_n_ss_2'.format(base_mod_dir)

def con_n_rch():
    n_con = mt.shape_file_to_model_array(
        "{}/inputs/N_inputs/shp/N_Conc_CP_adj_Rough.shp".format(sdp)
        , 'N_CMP_Conc')

    n_con_constant = deepcopy(n_con)
    n_con_constant[np.isfinite(n_con_constant)] = 10
    n_con_constant[np.isnan(n_con_constant)] = 0

    n_con[n_con > 200] = 200
    n_con[np.isnan(n_con)] = 0

    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': n_con_constant, 'ssm_stress_period_data': None,
        'adv_sov': 0, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 1, 'perlen': 1, 'nstp': 1, 'tsmult': 1,
        'ssflag': ['SState'],
        'dt0': 0, 'mxstrn': 10000, 'ttsmult': 1, 'ttsmax': 3.5,
        'gcg_isolve':3, 'gcg_inner':5000, 'gcg_outer':100
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/rch_constant_con'.format(base_dir),
            log = "C:/Users/matth/Desktop/constant_con_rch_ss.txt", safe_mode = False)

if __name__ == '__main__':
    #so far this is not working! The model may be too complex
    con_n_rch()