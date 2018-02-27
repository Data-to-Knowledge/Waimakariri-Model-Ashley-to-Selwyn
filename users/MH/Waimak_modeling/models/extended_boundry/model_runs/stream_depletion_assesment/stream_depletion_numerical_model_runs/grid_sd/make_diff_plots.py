# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 9:55 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import os
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt

def make_dif_plots(flux_8640_path, flux_432_path, outdir):
    flux_6840 = nc.Dataset(flux_8640_path)
    flux_432 = nc.Dataset(flux_432_path)

    for var in flux_6840.variables.keys():
        if var in ['longitude', 'latitude', 'layer', 'crs']:
            continue

        varoutdir = os.path.join(outdir, var)
        if not os.path.exists(varoutdir):
            os.makedirs(varoutdir)

        temp_6840 = np.array(flux_6840.variables[var])
        temp_432 = np.array(flux_432.variables[var])
        temp = temp_6840 - temp_432
        for layer in range(smt.layers - 1):
            vmin_opt = np.nanmin(temp[layer])
            vmax_opt = np.nanmax(temp[layer])
            bound = np.abs([vmin_opt, vmax_opt]).max()


            fig, ax = smt.plt_matrix(temp[layer], cmap='RdBu', base_map=True, vmin=bound*-1,vmax=bound,
                                     title='{} layer {} difference for sd (flux 6840 - flux 432)'.format(var, layer),
                                     no_flow_layer=layer)
            fig.savefig(os.path.join(varoutdir, 'layer_{:2d}_{}_difference.png'.format(layer, var)))
            plt.close(fig)

if __name__ == '__main__':
    path_432 = r"K:\mh_modeling\StrOpt_grid_sd\interpolated_StrOpt_sd_grid_data_flux_-432.0.nc"
    path_8640 = r"K:\mh_modeling\StrOpt_grid_sd\interpolated_StrOpt_sd_grid_data_flux_-8640.0.nc"
    temp_outdir = os.path.dirname(path_8640)
    model_id = os.path.basename(temp_outdir).split('_')[0]
    outdir = os.path.join(temp_outdir,'{}_8640_432_diff_plots'.format(model_id))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    make_dif_plots(path_8640,path_432,outdir)