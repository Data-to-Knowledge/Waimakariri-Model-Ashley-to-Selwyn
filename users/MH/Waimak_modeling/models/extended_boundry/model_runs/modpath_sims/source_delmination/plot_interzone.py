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

if __name__ == '__main__':
    #todo functionize and run for all of them
    data = nc.Dataset(
        r"T:\Temp\temp_gw_files\linwood_east_tla.nc")
    fig, axes = plt.subplots(2, 2, figsize=(18.5, 18.5))
    fig.suptitle('linwood_east_tla')
    model_xs, model_ys = smt.get_model_x_y()

    for i, name in enumerate(['forw_we_number', 'forw_st_number', 'back_st_number', 'back_we_number']):
        cur_map = np.array(data.variables[name])
        ax = axes.flatten()[i]
        smt.plt_matrix(smt.get_empty_model_grid()*np.nan, color_bar=False, base_map=True, ax=ax)
        levels = [0,1,2,5,10,15,20,30,40,50,75,100,165]
        norm = BoundaryNorm(levels, 256)
        c = ax.contourf(model_xs, model_ys, cur_map, alpha=0.5, levels=levels, cmap='RdYlBu', norm=norm)  # todo maybe spike alpha up
        ax.set_title(name)
    fig.colorbar(c, ticks=levels) #todo label all layers
    fig.savefig(r"C:\Users\MattH\Downloads\test_interzone.png") #todo check if pdf or png is smaller and easier to use

