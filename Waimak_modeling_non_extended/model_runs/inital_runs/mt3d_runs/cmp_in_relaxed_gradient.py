"""
Author: matth
Date Created: 31/05/2017 3:51 PM
"""

from __future__ import division
import numpy as np
import os
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
from users.MH.Waimak_modeling.supporting_data_path import sdp


if __name__ == '__main__':
    run_type = 2
    if run_type==2:
        model_type = 'm_relax_vert'
        n_con = mt.shape_file_to_model_array(
            "{}/inputs/N_inputs/shp/N_Conc_CP_adj_Rough.shp".format(sdp)
            , 'N_CMP_Conc')

        n_con[n_con > 200] = 200
        n_con[np.isnan(n_con)] = 0

        base_dir = "{}/mt3d_relaxed_vert".format(base_mod_dir)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        mt3d_input_cmp_con = {
            'm': None, 'mt3d_name': 'cmp',
            'ssm_crch': n_con, 'ssm_stress_period_data': None,
            'adv_sov': -1, 'adv_percel': 1,
            'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
            'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
            'nper': 10, 'perlen': 3650, 'nstp': [20, 20, 20, 20, 20, 20, 20, 20, 20, 20], 'tsmult': 1.1,
            'ssflag': None,
            'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
        }

        mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_cmp_con, '{}/{}_cmp_con'.format(base_dir, model_type),
                                       model_version=model_type)
