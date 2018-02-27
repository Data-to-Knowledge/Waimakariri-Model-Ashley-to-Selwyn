# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/02/2018 2:11 PM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt

if __name__ == '__main__':
    nc_data = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\mednload_unc.nc") #fix path

    for layer in range(smt.layers):
        data = np.array(nc_data.variables['mednload'][:,layer]) #this takes some memory ~ 13gb alternatively you can slice it by layers to arrive at the mean
    # netcdf slicing for the first layer and all realisations would be np.array(nc_data.variables['mednload'][:,0])  the dimensions (in the metadata) are (realisation, layer, row, col)
    # note it takes about the same amount of time to pull data out of a compressed netcdf regardless of the amount of data you pull out.
        mean = np.nanmean(data, axis=0)  # this is now a (k,i,j) array of the mean of all realisations n concentration
        std = np.nanstd(data, axis=0) # this is now a (k,i,j) array of the std of all realisations n concentration
        fig, ax = smt.plt_matrix(mean,title='Mean N Concentration Layer {}'.format(layer+1),no_flow_layer=layer,base_map=True)

    # my plotting mechanism is build into the infastructure of my scripts so is kind of tricky to package up, but I'm sure you've worked something up to quickly plot 2d arrays of this grid
