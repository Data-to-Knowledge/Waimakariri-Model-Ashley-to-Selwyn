# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 30/08/2017 11:25 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd

parameters = {  # note that these are in m3/s or %
    'pump_c': {'sd': 0.0682,  # fraction
               'min': np.nan,
               'max': np.nan,
               'inital': 1,
               'units': 'fraction'},

    'pump_s': {'sd': 0.0403,  # fraction
               'min': np.nan,
               'max': np.nan,
               'inital': 1,
               'units': 'fraction'},

    'pump_w': {'sd': 0.0547,  # fraction
               'min': np.nan,
               'max': np.nan,
               'inital': 1,
               'units': 'fraction'},

    'sriv': {'sd': 0.10,  # fraction
             'min': np.nan,
             'max': np.nan,
             'inital': 1,
             'units': 'fraction'},

    'n_race': {'sd': 0.25,  # fraction
               'min': np.nan,
               'max': np.nan,
               'inital': 1,
               'units': 'fraction'},

    'nbndf': {'sd': 0.255,  # fraction
              'min': np.nan,
              'max': np.nan,
              'inital': 1,
              'units': 'fraction'},

    's_race': {'sd': 0.2,  # fraction
               'min': np.nan,
               'max': np.nan,
               'inital': 1,
               'units': 'fraction'},

    'ulrzf': {'sd': 2/3,  # m3/s this is the total value it is distributed in my well adjust script
              'min': np.nan,
              'max': np.nan,
              'inital': 3,
              'units': 'm3/s'},

    'llrzf': {'sd': 1,  # m3/s this is the total value it is distributed in my well adjust script
              'min': np.nan,
              'max': np.nan,
              'inital': 0,
              'units': 'm3/s'},

    'top_e_flo': {'sd': 1.99 * 0.15,  # m3/s assumes 15% sd
                  'min': np.nan,
                  'max': np.nan,
                  'inital': 1.99,
                  'units': 'm3/s'},

    'top_c_flo': {'sd': 0.17 * 0.15,  # m3/s assumes 15% sd
                  'min': np.nan,
                  'max': np.nan,
                  'inital': 0.17,
                  'units': 'm3/s'},

    'mid_c_flo': {'sd': 0.04 * 0.15,  # m3/s assumes 15% sd
                  'min': np.nan,
                  'max': np.nan,
                  'inital': 0.04,
                  'units': 'm3/s'},

    'rch_wai_ir': {'sd': 71 / 308,  # fraction
                   'min': 94 / 308,
                   'max': 521 / 308,
                   'inital': 1,
                   'units': 'fraction'},

    'rch_sel_ir': {'sd': 54 / 250,  # fraction
                   'min': 88 / 250,
                   'max': 412 / 250,
                   'inital': 1,
                   'units': 'fraction'},

    'rch_drylnd': {'sd': 44 / 192,  # fraction
                   'min': 61 / 192,
                   'max': 301 / 192,
                   'inital': 1,
                   'units': 'fraction'},

    'rch_conf': {'sd': 20 / 94,  # fraction
                 'min': 32 / 94,
                 'max': 154 / 94,
                 'inital': 1,
                 'units': 'fraction'}
}


def get_param_table(convert_to_m3day = True):
    outdata = pd.DataFrame(parameters).transpose()
    idx = pd.isnull(outdata.loc[:, 'max'])
    outdata.loc[idx, 'max'] = 3 * outdata.loc[idx, 'sd'] + outdata.loc[idx, 'inital']
    idx = pd.isnull(outdata.loc[:,'min'])
    outdata.loc[idx, 'min'] = -3 * outdata.loc[idx, 'sd'] + outdata.loc[idx, 'inital']

    if convert_to_m3day:
        idx = outdata.units=='m3/s'
        outdata.loc[idx,[u'inital', u'max', u'min', u'sd']] *= 86400
        outdata.loc[idx,'units'] = 'm3/day'
    return outdata


def create_parameter_file(unc_file_path):
    # inital lines just a place holder for now updated by cath
    inital_lines = [
        'START COVARIANCE_MATRIX\n',
        'file hk_cov.mat\n',
        'variance_multiplier 5.0\n',
        'END COVARIANCE_MATRIX\n',
        '\n',
        'START COVARIANCE_MATRIX\n',
        'file por_cov.mat\n',
        'variance_multiplier 5.0\n',
        'END COVARIANCE_MATRIX\n',
        '\n'
    ]

    normal_actual_parameters = ['top_e_flo', 'top_c_flo', 'mid_c_flo','llrzf', 'ulrzf']
    normal_mult_parameters = ['pump_c', 'pump_s', 'pump_w', 'sriv', 'n_race', 'nbndf','s_race']

    uniform_actual_parameters = []  # no uniforms used
    uniform_mult_parameters = []

    with open(unc_file_path, 'w') as f:
        f.writelines(inital_lines)
        # write normal actual parameters
        f.writelines([
            'START STANDARD_DEVIATION\n',
            'std_multiplier 1\n'
        ])

        # write normal actual paramters
        for param in normal_actual_parameters:
            sd = parameters[param]['sd']
            _min = parameters[param]['min']
            _max = parameters[param]['max']
            if np.isnan(sd) or np.isfinite(_min) or np.isfinite(_max):
                raise ValueError('{} has a nan or actual value where it should not'.format(param))
            f.write('{} {}\n'.format(param, sd * 86400))

        # write normal multiplicitve paramters
        for param in normal_mult_parameters:
            sd = parameters[param]['sd']
            _min = parameters[param]['min']
            _max = parameters[param]['max']
            if np.isnan(sd) or np.isfinite(_min) or np.isfinite(_max):
                raise ValueError('{} has a nan or actual value where it should not'.format(param))
            f.write('{} {}\n'.format(param, sd))

        f.write('END STANDARD_DEVIATION\n')
        f.write('\n')

        # write uniform actual parameters never used this
        for param in uniform_actual_parameters:
            sd = parameters[param]['sd']
            _min = parameters[param]['min']
            _max = parameters[param]['max']
            if np.isfinite(sd) or np.isnan(_min) or np.isnan(_max):
                raise ValueError('{} has a nan or actual value where it should not'.format(param))
            f.write('{} {} {}\n'.format(param, _min * 86400, _max * 86400))

        for param in uniform_mult_parameters:
            sd = parameters[param]['sd']
            _min = parameters[param]['min']
            _max = parameters[param]['max']
            if np.isfinite(sd) or np.isnan(_min) or np.isnan(_max):
                raise ValueError('{} has a nan or actual value where it should not'.format(param))
            f.write('{} {} {}\n'.format(param, _min, _max))


if __name__ == '__main__':
    test = get_param_table()
    test2 = get_param_table(False)
    print'done'