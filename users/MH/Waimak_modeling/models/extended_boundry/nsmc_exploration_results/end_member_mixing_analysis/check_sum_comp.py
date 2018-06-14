# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/11/2017 7:09 PM
"""

from __future__ import division
from core import env
import flopy_mh as flopy
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt

alpine = flopy.utils.UcnFile(env.gw_met_data(r"mh_modeling\data_from_gns\mt_aw_ex_alpine_philow\mt_aw_ex_alpine_philow.ucn")).get_data(kstpkper=(0,0))
inland = flopy.utils.UcnFile(env.gw_met_data(r"mh_modeling\data_from_gns\EM_inland_ucnrepo\mt_aw_ex_inland_philow.ucn")).get_data(kstpkper=(0,0))
coastal = flopy.utils.UcnFile(env.gw_met_data(r"mh_modeling\data_from_gns\EM_coast_ucnrepo\mt_aw_ex_coastal_philow.ucn")).get_data(kstpkper=(0,0))

alpine[alpine==-1] = np.nan
inland[inland==-1] = np.nan
coastal[coastal==-1] = np.nan

sums = alpine+inland+coastal
print ('min: {}'.format(np.nanmin(sums)))
print ('1st: {}'.format(np.nanpercentile(sums,1)))
print ('5th: {}'.format(np.nanpercentile(sums,5)))
print ('25th: {}'.format(np.nanpercentile(sums,25)))
print ('50th: {}'.format(np.nanpercentile(sums,50)))
print ('75th: {}'.format(np.nanpercentile(sums,75)))
print ('95th: {}'.format(np.nanpercentile(sums,95)))
print ('99th: {}'.format(np.nanpercentile(sums,99)))
print ('max: {}'.format(np.nanmax(sums)))

for l in range(smt.layers):
    smt.plt_matrix(sums[l],title='layer {}'.format(l))

plt.show()