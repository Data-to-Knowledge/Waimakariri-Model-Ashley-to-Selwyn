# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/12/2017 12:43 PM
"""

from __future__ import division

import os
import pickle
from time import time

import numpy as np

from con_array import _make_mednload_approx
from run_emulator import run_emulator
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.extract_data import save_forward_data
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.setup_forward_modpath import setup_run_forward_modpath
from waimak_extended_boundry.model_run_tools import get_cbc

if __name__ == '__main__':
    #just a test
    run_model = True
    t = time()
    mp_ws = r"D:\mh_waimak_models\modpath_emulator\try_truncation"
    mp_name = 'NsmcBase_first_try'
    path_file = os.path.join(mp_ws, 'NsmcBase_first_try.mppth')
    if run_model:
        if not os.path.exists(mp_ws):
            os.makedirs(mp_ws)
        cbc = get_cbc('NsmcBase', mp_ws)
        print('{} min to make cbc'.format((time() - t) / 60))
        t = time()

        setup_run_forward_modpath(cbc, mp_ws, mp_name, min_part=1, max_part=100)
        print('{} min to setup and run modpath'.format((time() - t) / 60))
        t = time()

        save_forward_data(path_file, path_file.replace('.mppth', '.hdf'))
        print('{} min to make emulator'.format((time() - t) / 60))
        t = time()

    bnd_type = np.loadtxt(os.path.join(mp_ws,'{}_bnd_type.txt'.format(mp_name)))
    load = _make_mednload_approx(bnd_type)
    outdata = run_emulator(path_file.replace('.mppth','.hdf'),load, bd_type=bnd_type)
    pickle.dump(outdata, open(os.path.join(os.path.dirname(path_file),'{}_test_med_n.p'.format(mp_name)),'w'))
    print('{} min to run emulator'.format((time() - t) / 60))


