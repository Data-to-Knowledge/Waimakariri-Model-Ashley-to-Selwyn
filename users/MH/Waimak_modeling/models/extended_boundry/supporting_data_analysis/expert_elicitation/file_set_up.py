# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/08/2017 9:10 AM
"""

from __future__ import division
from core import env
import os
import pandas as pd
from copy import deepcopy

users = [
         #'paul_white',
         'peter_callander',
         'scott_wilson',
         'jens_rekker',
         'john_talbot',
         #'lee_burbery',
         'zeb_etheridge']

user_codes = {'jens_rekker': 'A',
              'john_talbot': 'B',
              'lee_burbery': 'C',
              'peter_callander': 'D',
              'paul_white': 'E',
              'scott_wilson': 'F',
              'zeb_etheridge': 'G'}

major_qoi = ['LSR', 'pumping', 'LRZ_flux', 'offshore', 'possible_ks', 'race_loss']
minor_qoi = {'pumping': ['model_average_pumping'], #, 'extra_1', 'extra_2', 'extra_3'],
             'possible_ks': ['possible_ks'], #'extra_1', 'extra_2', 'extra_3'],
             'offshore': ['N_waimak', 'S_waimak', 'CHCH', 'Selwyn'], #  'extra_1', 'extra_2', 'extra_3'],
             'LRZ_flux': ['LRZ_flux'],# 'extra_1', 'extra_2', 'extra_3'],
             'LSR': ['waimak_irrigated', 'waimak_dryland', 'inland_chch',
                     'coastal_chch', 'selwyn_irrigated','selwyn_dryland'],#,
                     #'extra_1', 'extra_2', 'extra_3'],
             'race_loss': ['waimak_races','sewlyn_races']# 'extra_1', 'extra_2', 'extra_3']
            }
minor_qoi_shapes = {'pumping': (1,1),
             'possible_ks': (1,1),
             'offshore': (2,2),
             'LRZ_flux': (1,1),
             'LSR': (2,3),
             'race_loss': (2,1)}

values = ['lower_limit', 'upper_limit', 'median', 'lower_quartile', 'upper_quartile']

def make_empty_files(base_dir,date, include_extras=True):
    if os.path.exists(base_dir):
        raise ValueError('the path already exists, do not overwrite!')
    os.makedirs(base_dir)
    for u in users:
        os.makedirs('{}/{}'.format(base_dir,u))
        for qoi in major_qoi:
            index = deepcopy(minor_qoi[qoi])
            if include_extras:
                index.extend(['extra_1', 'extra_2', 'extra_3'])
            temp = pd.DataFrame(columns=values,index=index)
            temp.to_excel('{b}/{us}/W_elicitation_{d}_{us}_{q}.xlsx'.format(b=base_dir, us=u, d=date, q=qoi))



if __name__ == '__main__':
    make_empty_files(r"C:\Users\MattH\Desktop\test_waimak_elicitation5",'21-08-2017')

    print 'done'