# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/03/2018 3:15 PM
"""

from __future__ import division
from core import env
from hydraulics_cleaned import hunt2003
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import itertools


# is it worth doing a bulk comparison for each stream (e.g. all reaches?)
# I should consider dry reaches as well
# point to point nodes for all shapefiles might be the simplest
# consider an assumed pumped aquifer thickness of 25m starting from the well midpoint.  area above is the 'aquitard'

#todo resort out how to handle indexes (maybe vectorize later...)

def _mean_of_items(idxs, array):
    np.nanmean(array[idxs])

def vec_mean_of_items(indexes, array):
    return np.vectorize(_mean_of_items,excluded=['array'])(indexes,array)

def get_indexes_2d(row, col, radius):  # note when vectorizing return set otyps to object
    num_cells = radius // smt.grid_space
    rows = np.arange(row - num_cells, row + num_cells + 1)
    rows = rows[(rows >= 0) & (rows < smt.rows)]  # get rid of all rows, and cols which are out of order
    cols = np.arange(col - num_cells, col + num_cells + 1)
    cols = cols[(cols >= 0) & (cols < smt.cols)]
    idxs = np.array(list(itertools.product(rows, cols)))
    idxs = [idxs[:, 0], idxs[:, 1]]
    return idxs


def get_tranmissivity(indexes_2d, layer, all_kh):
    """
    get the transmissivity as an average of teh area
    :param indexes_2d: all of the 2-d indexes
    :param layer: well layer
    :param radius: radius to consider (m) this is considered the square rounded number of cells that the radius
                   encounters
    :return:
    """

    # get ks
    ks = vec_mean_of_items(indexes_2d, all_kh[layer])

    return ks*25


def get_aquitard_kv(indexes_2d, pumping_layer, all_kv):
    """
    get the kv of the aquitard as an average of teh area
    :return:
    """
    layer_kvs = []
    for layer in range(pumping_layer):
        layer_kvs.append(vec_mean_of_items(indexes_2d,all_kv[layer]))
    return np.nanmean(layer_kvs)

def get_aquitard_sat_thickness(indexes_2d, pumping_layer):
    """
    """
    # what do i need to do here for layer 0
    raise NotImplementedError


def get_separation_distance():
    raise NotImplementedError
