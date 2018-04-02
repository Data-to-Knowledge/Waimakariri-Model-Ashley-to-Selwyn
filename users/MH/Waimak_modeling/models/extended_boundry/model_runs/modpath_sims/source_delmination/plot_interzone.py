# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 28/03/2018 8:00 PM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import os
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt


def plot_interzone_nc(path):
    data = nc.Dataset(path)
    fig, axes = plt.subplots(2, 2, figsize=(18.5, 18.5))
    fig.suptitle(os.path.basename(path).replace('.nc', ''))
    model_xs, model_ys = smt.get_model_x_y()
    levels = [0, 1, 2, 5, 10, 15, 20, 30, 40, 50, 75, 100, 165]
    norm = BoundaryNorm(levels, 256)

    for i, name in enumerate(['forw_st_number', 'back_st_number', 'forw_we_number', 'back_we_number']):
        cur_map = np.array(data.variables[name])
        ax = axes.flatten()[i]
        smt.plt_matrix(smt.get_empty_model_grid() * np.nan, color_bar=False, base_map=True, ax=ax)
        c = ax.contourf(model_xs, model_ys, cur_map, alpha=0.5, levels=levels, cmap='RdYlBu',
                        norm=norm)
        ax.set_title(name)
    fig.colorbar(c, ticks=levels)
    return fig


if __name__ == '__main__':
    for rt in ['stocastic set', 'AshOpt']:
        base_outdir = env.gw_met_data(r"mh_modeling\interzone_source_zones\{}\plots".format(rt))
        if not os.path.exists(base_outdir):
            os.makedirs(base_outdir)
        paths = glob(env.gw_met_data(r"mh_modeling\interzone_source_zones\{}\individual_netcdfs\*.nc".format(rt)))
        for path in paths:
            fig = plot_interzone_nc(path)
            fig.savefig(os.path.join(base_outdir,os.path.basename(path).replace('.nc', '.png')))
