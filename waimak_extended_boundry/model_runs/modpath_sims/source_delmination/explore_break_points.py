# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/02/2018 8:35 AM
"""

from __future__ import division
import numpy as np
import netCDF4 as nc

if __name__ == '__main__':

    remove_keys = {'crs', 'latitude', 'longitude'}

    # single sites
    backs = []
    # backward strong
    bs = nc.Dataset(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\single_source_zones\stocastic set\backward_strong.nc")
    variables = list(set(bs.variables.keys())-remove_keys)
    for key in variables:
        temp = np.array(bs.variables[key])
        temp2 = temp[~np.isclose(temp, 0)]
        backs.append(temp2)
    backs = np.concatenate(backs).flatten()

    fors = []
    # forward strong
    fs = nc.Dataset(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\single_source_zones\stocastic set\backward_strong.nc")
    variables = list(set(fs.variables.keys())-remove_keys)
    for key in variables:
        temp = np.array(fs.variables[key])
        temp2 = temp[~np.isclose(temp, 0)]
        fors.append(temp2)
    fors = np.concatenate(fors).flatten()


    print('done') # I am playing with histograms no great thoughts ... I think it is a qualitative, subjective thing.


