# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 28/09/2017 5:02 PM
"""

from __future__ import division
from core import env
import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt


if __name__ == '__main__':
    # a test script to look at teh irrigation demand
    lsrm_rch_base_dir = env.gw_met_data('niwa_netcdf/lsrm/lsrm_results/water_year_means')
    paths = glob(os.path.join(lsrm_rch_base_dir,'*/ird*.txt'))

    data = []
    for path in paths:
        with open(path) as f:
            data.append(float(f.readline()))
    data = np.array(data)

    print 'done'
