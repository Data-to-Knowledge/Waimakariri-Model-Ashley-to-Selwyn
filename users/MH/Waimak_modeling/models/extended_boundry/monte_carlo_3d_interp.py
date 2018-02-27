# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/07/2017 2:31 PM
"""

from __future__ import division
from core import env
from scipy.interpolate import Rbf
import pandas as pd
import multiprocessing
import logging
import numpy as np
from copy import deepcopy

#this did not work well, it's unclear if it was just our data or the method
def interp_3d_montecarlo(locations, data, method, num_sim=None, return_values=False):
    """
    returns values from a montecarloed interpolation from the data provided (uses radial basis function interpolation)
    pulls montecarlos from a normal distribution defined as mean=d, and std = sd
    :param locations: data_frame with(x,y,z)
    :param data: dataframe with x,y,z,d,sd,
    :param method: which interpolation method to use
    :param num_sim: int the number of monte-calro simulations
    :return_values: bool if True also returns a list of values for each point (as a dataframe column)
    :return: dataframe (x,y,z,median,sd,[values])
    """
    outdata = deepcopy(locations)

    # define runs
    runs=[None]*num_sim

    for i in range(num_sim):
        runs[i]={'x':data['x'],
          'y':data['y'],
          'z': data['z'],
          'd':np.random.normal(data['d'],data['sd']),
          'locations': (np.array(locations['x']),np.array(locations['y']),np.array(locations['z'])),
          'method': method}


    pool_size = multiprocessing.cpu_count()-1
    pool = multiprocessing.Pool(processes=pool_size)
    pool_outputs = pool.map(_interp_3d_mp, runs)
    pool.close()  # no more tasks
    pool.join()

    outdata['data_mean'] = np.array(pool_outputs).mean(axis=0)
    outdata['data_sd'] = np.array(pool_outputs).std(axis=0)

    if return_values:
        values = pd.DataFrame(np.array(pool_outputs).transpose())
        values.index.names = ['site_index']
        return outdata, values
    else:
        return outdata


def _interp_3d_mp(kwargs):
    return interp_3d(**kwargs)

def interp_3d(x,y,z,d,locations, method):
    """
    return values from a 3d interpolation (uses radial basis function interpolation)

    :param x: array
    :param y: array
    :param z: array
    :param d: array
    :param locations: tuple of (x,y,z) arrays
    :param method: interpolation method to use
    :return:
    """

    interp = Rbf(x,y,z,d, function=method)
    return interp(*locations)


