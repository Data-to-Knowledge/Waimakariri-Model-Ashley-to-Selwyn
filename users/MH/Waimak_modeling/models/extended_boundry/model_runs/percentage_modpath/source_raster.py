# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/12/2017 9:56 AM
"""

from __future__ import division
from core import env
from time import time
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

def define_source(emulator_path, bd_type, index): # todo could make this into a list of indexes or dictionary of indexes...
    """

    :param emulator_path: path to the emulator (hdf)
    :param index: a boolean array of size smt.layers, rows, cols which represent the area of interest
    :return: array of shape (smt.layers, rows, cols) with a percentage of water from source
    """
    # run some checks on inputs
    if index is None:
        index = smt.get_no_flow() == 1 # all active cells
    elif isinstance(index, np.ndarray):
        if index.shape != (smt.layers,smt.rows,smt.cols) or index.dtype != bool:
            raise ValueError('index must have shape {} and be boolean'.format((smt.layers,smt.rows,smt.cols)))
    else:
        raise ValueError('index must be None or ndarray not {}'.format(type(index)))

    # load emulator and initialize outdata
    print('loading emulator')
    t = time()
    emulator = pd.read_hdf(emulator_path)  # this keeps the structure of everything

    outdata = smt.get_empty_model_grid(False)
    print('took {} s to load emulator'.format(time()-t))
    t = time()


    # get area of interest
    print('calculating area of interest')
    layers, rows, cols = np.meshgrid(range(smt.layers),range(smt.rows), range(smt.cols),indexing='ij')
    ids = ['{:02d}_{:03d}_{:03d}'.format(k, i, j) for k, i, j in zip(layers[index],rows[index],cols[index])]
    temp = np.in1d(emulator.index.values,ids)
    emulator = emulator.loc[temp]
    print('took {} s to identify area'.format(time()-t))
    t = time()

    # calculate source percentage
    emulator = emulator.reset_index()
    temp = emulator.groupby('Particle_Group').aggregate({'fraction':np.sum})
    temp *= 1/temp.sum()

    # populate array
    ibnd = smt.get_no_flow(0).flatten()
    idx = bd_type.flatten() != -1
    outdata = outdata.flatten()
    temp_array = outdata[idx]
    temp_array[temp.index.values-1] = temp.fraction.values
    outdata[idx] = temp_array
    outdata = outdata.reshape((smt.rows,smt.cols))
    smt.plt_matrix(outdata!=0)
    range(1, idx.sum() + 1)

    #todo something is wrong with indexing....  work bd type through
    return outdata

if __name__ == '__main__':
    from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
    index = smt.get_empty_model_grid(True)
    index[:,150:171,290:301] = True
    define_source(r"T:\Temp\temp_gw_files\first_try.hdf",index.astype(bool))