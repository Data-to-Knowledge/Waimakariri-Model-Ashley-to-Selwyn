"""
Author: matth
Date Created: 13/05/2017 8:21 AM
"""

from __future__ import division
from future.builtins import input
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
import pandas as pd
import os
import flopy
import logging
import numpy as np
from copy import deepcopy
import multiprocessing
from users.MH.Waimak_modeling.supporting_data_path import sdp

base_dir = '{}/component_con_n_zero_start'.format(base_mod_dir) # add to notes
if not os.path.exists(base_dir):
    os.makedirs(base_dir)
itype = flopy.mt3d.Mt3dSsm.itype_dict()
# constant N recharge

def con_n_rch_N():
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
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/rch_constant_con'.format(base_dir),
            log = "C:/Users/matth/Desktop/constant_con_rch2.txt", safe_mode = False)

def con_n_rch_S():
    n_con = mt.shape_file_to_model_array(
        "{}/inputs/N_inputs/shp/N_Conc_CP_adj_Rough.shp".format(sdp)
        , 'N_CMP_Conc')

    n_con_s = mt.shape_file_to_model_array("{}/inputs/shp_files/rch_s/rch_s.shp".format(sdp),
                                           'nload')

    n_con_s[np.isfinite(n_con)] = 0
    n_con_s[np.isnan(n_con_s)] = 0
    n_con_s[0:94,300:] = 0

    n_con_constant = n_con_s


    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': n_con_constant, 'ssm_stress_period_data': None,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/rch_S_constant_con'.format(base_dir),
            log = "C:/Users/matth/Desktop/constant_con_rch2_S.txt", safe_mode = False)

# constant N for streams (one at a time? yes)

def con_n_wai():
    segs = list(set(mt.get_all_segs(35)) - set(mt.get_all_segs(17)))
    str_data = mt.wraps.mflow.get_base_mf_ss().str.stress_period_data.data[0]
    str_data = str_data[np.in1d(str_data['segment'],segs)]
    stress_data = pd.DataFrame()
    stress_data['k'] = str_data['k']
    stress_data['i'] = str_data['i']
    stress_data['j'] = str_data['j']
    stress_data['css'] = 10
    stress_data['itype'] = 21
    stress_data = {0: list(stress_data.values)} # this might cause problems later with different data types ssee if can pass DF

    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': np.zeros((190,365)), 'ssm_stress_period_data': stress_data,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/wai_constant_con'.format(base_dir),
                                   log="C:/Users/matth/Desktop/constant_con_wai2.txt", safe_mode=False)


def con_n_ashley():
    segs = mt.get_all_segs(34)
    str_data = mt.wraps.mflow.get_base_mf_ss().str.stress_period_data.data[0]
    str_data = str_data[np.in1d(str_data['segment'],segs)]
    stress_data = pd.DataFrame()
    stress_data['k'] = str_data['k']
    stress_data['i'] = str_data['i']
    stress_data['j'] = str_data['j']
    stress_data['css'] = 10
    stress_data['itype'] = 21
    stress_data = {0: list(stress_data.values)}
    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': np.zeros((190,365)), 'ssm_stress_period_data': stress_data,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/ash_constant_con'.format(base_dir),
                                   log="C:/Users/matth/Desktop/constant_con_ashley2.txt", safe_mode=False)


def con_n_cust():
    segs = mt.get_all_segs(37)
    str_data = mt.wraps.mflow.get_base_mf_ss().str.stress_period_data.data[0]
    str_data = str_data[np.in1d(str_data['segment'],segs)]
    stress_data = pd.DataFrame()
    stress_data['k'] = str_data['k']
    stress_data['i'] = str_data['i']
    stress_data['j'] = str_data['j']
    stress_data['css'] = 10
    stress_data['itype'] = 21
    stress_data = {0: list(stress_data.values)}
    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': np.zeros((190,365)), 'ssm_stress_period_data': stress_data,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/cust_constant_con'.format(base_dir),
                                   log="C:/Users/matth/Desktop/constant_con_cust2.txt", safe_mode=False)


def con_n_eyre():
    segs = mt.get_all_segs(17)
    str_data = mt.wraps.mflow.get_base_mf_ss().str.stress_period_data.data[0]
    str_data = str_data[np.in1d(str_data['segment'],segs)]
    stress_data = pd.DataFrame()
    stress_data['k'] = str_data['k']
    stress_data['i'] = str_data['i']
    stress_data['j'] = str_data['j']
    stress_data['css'] = 10
    stress_data['itype'] = 21
    stress_data = {0: list(stress_data.values)}

    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': np.zeros((190,365)), 'ssm_stress_period_data': stress_data,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/eyre_constant_con'.format(base_dir),
                                   log="C:/Users/matth/Desktop/constant_con_eyre2.txt", safe_mode=False)

#constant N for injection wells
# debug
def con_n_races():
    well_data = mt.wraps.mflow.get_base_mf_ss().wel.stress_period_data.data[0]
    well_data = pd.DataFrame(well_data[well_data['flux'] > 0])
    well_data = well_data.drop('iface',axis=1)
    well_data['flux'] = 10 # set a constant concentration
    well_data['itype'] = itype['WEL']

    stress_data = {0: list(well_data.values)}
    mt3d_input_constant_con = {
        'm': None, 'mt3d_name': 'constant_concentration',
        'ssm_crch': np.zeros((190,365)), 'ssm_stress_period_data': stress_data,
        'adv_sov': -1, 'adv_percel': 1,
        'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
        'dsp_lon': 0.0005, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,  # this line set
        'nper': 10, 'perlen': 3650, 'nstp': 20, 'tsmult': 1.1,
        'ssflag': None,
        'dt0': [0.25, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], 'mxstrn': 10000, 'ttsmult': 1.1, 'ttsmax': 180
    }

    mt.wraps.mt3d.run_modflow_mt3d(mt3d_input_constant_con, '{}/races_constant_con'.format(base_dir),
                                   log="C:/Users/matth/Desktop/constant_con_races2.txt", safe_mode=False)


if __name__ == '__main__':
    rtype = 1
    if rtype == 1:
        continue_ = input("continue may overwrite many files, for constant concentration runs y/n")
        if continue_.lower() != 'y':
            raise ValueError('aborted to prevent overwrite')

        jobs = []
        multiprocessing.log_to_stderr(logging.DEBUG)
        #m_run = multiprocessing.Process(target=con_n_rch_N, name='rch')
        #jobs.append(m_run)
        #m_run.start()

        m_run2 = multiprocessing.Process(target=con_n_wai, name='wai')
        jobs.append(m_run2)
        m_run2.start()

        m_run3 = multiprocessing.Process(target=con_n_ashley, name='ash')
        jobs.append(m_run3)
        m_run3.start()

        m_run4 = multiprocessing.Process(target=con_n_cust, name='cust')
        jobs.append(m_run4)
        m_run4.start()

        m_run5 = multiprocessing.Process(target=con_n_eyre, name='eye')
        jobs.append(m_run5)
        m_run5.start()

        #m_run6 = multiprocessing.Process(target=con_n_races, name='race')
        #jobs.append(m_run6)
        #m_run6.start()

        m_run7 = multiprocessing.Process(target=con_n_rch_S(), name='RCH_S')
        jobs.append(m_run7)
        m_run7.start()

    elif rtype == 2:
        con_n_wai()


