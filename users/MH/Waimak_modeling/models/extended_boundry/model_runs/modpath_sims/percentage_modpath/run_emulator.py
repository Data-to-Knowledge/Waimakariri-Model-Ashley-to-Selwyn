# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/12/2017 6:10 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import numpy as np
import pandas as pd
from time import time

# make a function to run the full emulator and to run only certain cells given a load layer.  also function to run as a
# stocastic simulation

def get_group_number(index, bd_type):
    """

    :param index:
    :param bd_type:
    :return:
    """
    assert isinstance(index,np.ndarray)
    assert index.dtype == bool
    if index.shape != (smt.rows,smt.cols):
        raise ValueError('array must have shape of {} not {}'.format((smt.rows,smt.cols),index.shape))
    idx = bd_type.flatten() != -1
    temp = np.where(index.flatten()[idx])[0]
    return temp



def _convert_data_to_cell_dict(array, bd_type):
    """
    converts the data to a dictionary to replace teh group ids
    :param array: array of data (will be converted to float)
    :return: dict(cell_id: value)
    """
    assert isinstance(array,np.ndarray)
    if array.shape != (smt.rows,smt.cols):
        raise ValueError('array must have shape of {} not {}'.format((smt.rows,smt.cols),array.shape))
    array = array.astype(float)
    idx = bd_type.flatten() != -1
    outdata = {k:v for k,v in zip(np.array(range(1,idx.sum()+1))*-1, array.flatten()[idx])}
    return outdata


def run_emulator(emulator_path, load_layer, bd_type, index=None):
    """

    :param emulator_path: path to the emulator (hdf)
    :param load_layer:
    :param bd_type: the boundary condition type used to map particles
    :param index: a boolean array of size smt.layers, rows, cols or None,
                  if None concentrations will be calculated across all active cells
    :return:
    """
    # run some checks on inputs
    if index is None:
        index = smt.get_no_flow() == 1 # all active cells
    elif isinstance(index, np.ndarray):
        if index.shape != (smt.layers,smt.rows,smt.cols) or index.dtype != bool:
            raise ValueError('index must have shape {} and be boolean'.format((smt.layers,smt.rows,smt.cols)))
    else:
        raise ValueError('index must be None or ndarray not {}'.format(type(index)))
    assert isinstance(load_layer, np.ndarray)
    assert load_layer.shape == (smt.rows, smt.cols)

    # load emulator and initialize outdata
    print('loading emulator')
    t = time()
    emulator = pd.read_hdf(emulator_path)  # this keeps the structure of everything

    outdata = smt.get_empty_model_grid(True)
    outdata.fill(np.nan)
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

    # add concentrations
    print('calculating concentrations')
    cons = _convert_data_to_cell_dict(load_layer, bd_type)
    print('replacing cons')
    emulator['con'] = vec_translate(emulator.loc[:, 'Particle_Group'].values*-1, cons) # much faster and more effcient than replace
    print('sorting out data')
    emulator['con'] *= emulator['fraction']
    temp_out = emulator['con'].groupby('ref_cell_id').sum()

    print('took {} s to calculate concentrations'.format(time()-t))
    t = time()

    # output to an array
    outdata[index] = temp_out.loc[ids].values
    print('took {} s to output concentrations'.format(time()-t))

    # could either gausian blur or run a spline to to account for dispersion
    return outdata

# quickest way to assign is to loop or to use below
def vec_translate(a, d):
    return np.vectorize(d.__getitem__)(a)

