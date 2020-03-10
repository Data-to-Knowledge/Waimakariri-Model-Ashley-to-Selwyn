"""
Author: matth
Date Created: 11/05/2017 9:33 AM
"""

from __future__ import division
from future.builtins import input
import numpy as np
import os
import logging
import multiprocessing
from copy import deepcopy
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, sdp


def main():
    # load N layer
    n_con = mt.shape_file_to_model_array(
        "{}/inputs/N_inputs/shp/N_Conc_CP_adj_Rough.shp".format(sdp)
        , 'N_CMP_Conc')

    n_con_constant = deepcopy(n_con)
    n_con_constant[np.isfinite(n_con_constant)] = 10
    n_con_constant[np.isnan(n_con_constant)] = 0

    n_con[n_con > 200] = 200
    n_con[np.isnan(n_con)] = 0

    # set up mt3d options as dict
    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': n_con_constant, 'ssm_stress_period_data': None,
        'adv_sov': 0, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 1, 'btn_nprs': -10, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01, # this line set
        'nper': 10, 'perlen': 3650, 'nstp': [100,75,50,25,20,20,20,20,20,20], 'tsmult': 1.1, 'ssflag': None,
        'dt0': [0.25,0.3,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5], 'mxstrn': 1000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt3d_input_layer = {
        'm': None, 'mt3d_name': 'CMP',
        'ssm_crch': n_con, 'ssm_stress_period_data': None,
        'adv_sov': 0, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 1, 'btn_nprs': -10, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01, # this line set
        'nper': 10, 'perlen': 3650, 'nstp': [100,75,50,25,20,20,20,20,20,20], 'tsmult': 1.1, 'ssflag': None,
        'dt0': [0.25,0.3,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5], 'mxstrn': 1000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    #run both a constant load and the cell load in parallel
    base_dir = "{}/first_mt3d2".format(base_mod_dir)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    mt3d_inputs = [mt3d_input_constant_con, mt3d_input_layer]
    process_names = ['constant_con','layer']
    logs = ["C:/Users/matth/Desktop/constant_con_buff.txt","C:/Users/matth/Desktop/layer_buff.txt"]
    dirs = ['{}/constant_con'.format(base_dir),'{}/cmp_layer'.format(base_dir)]
    for dir_path in dirs:
        if os.path.exists(dir_path):
            cont = input(
                'run_modflow_mt3d will delete the directory:\n {} \n continue y/n\n'.format(dir_path)).lower()
            if cont == 'n':
                raise ValueError('script aborted so as not to overwrite {}'.format(dir_path))

    jobs = []
    for i in range(len(mt3d_inputs)):
        multiprocessing.log_to_stderr(logging.DEBUG)
        m_run = multiprocessing.Process(target=mt.wraps.mt3d.run_modflow_mt3d,
                                        args=(mt3d_inputs[i],dirs[i],logs[i],False),
                                        name=process_names[i])
        jobs.append(m_run)
        m_run.start()

if __name__ == '__main__':
    type = 1
    if type == 1:
        main()
    elif type==2:
        n_con = mt.shape_file_to_model_array(
            "{}/inputs/N_inputs/shp/N_Conc_CP_adj_Rough.shp".format(sdp)
            , 'N_CMP_Conc')

        n_con_constant = deepcopy(n_con)
        n_con_constant[np.isfinite(n_con_constant)] = 10
        n_con_constant[np.isnan(n_con_constant)] = 0

        n_con[n_con > 200] = 200
        n_con[np.isnan(n_con)] = 0

        base_dir = "{}/first_mt3d2".format(base_mod_dir)
        mt3d_input_constant_con = {
            'm': None, 'mt3d_name': 'constant_concentration',
            'ssm_crch': n_con_constant, 'ssm_stress_period_data': None,
            'adv_sov': 0, 'adv_percel': 1,
            'btn_porsty': 0.05, 'btn_scon': 1, 'btn_nprs': 0, 'btn_timprs': None,
            'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
            'nper': 10, 'perlen': 3650, 'nstp': [100, 75, 50, 25, 20, 20, 20, 20, 20, 20], 'tsmult': 1.1,
            'ssflag': None,
            'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 1000, 'ttsmult': 1.1, 'ttsmax': 180
        }

        mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con,'{}/cmp_constant_con'.format(base_dir))