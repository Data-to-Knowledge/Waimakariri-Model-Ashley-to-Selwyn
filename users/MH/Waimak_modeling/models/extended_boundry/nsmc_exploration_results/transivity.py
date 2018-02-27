# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 24/11/2017 12:54 PM
"""

from __future__ import division
from core import env
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import netCDF4 as nc
from nsmc_visualisation.parameters.spatially_average_param import plot_sd_mean_multid
import matplotlib.pyplot as plt

def _log10_t(x,layer,**kwargs):
    elv_db = smt.calc_elv_db()
    x = 10**x
    thickness = elv_db[0:-1] - elv_db[1:]
    outdata = np.log10(x*thickness[layer])

    return outdata

if __name__ == '__main__':
    nc_parm = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    fig, ax = plot_sd_mean_multid(filter_strs=['emma_no_wt', 'run_mt3d'], layer=10, nc_param_data=nc_parm, data_id='kh',
                        title='T for layer 11', basemap=True, contour={'sd': False, 'mean': False},
                        contour_color='g', vmins=None, vmaxes=None, tranform=_log10_t)

    plt.show()