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

# is it worth doing a bulk comparison for each stream (e.g. all reaches?)
# I should consider dry reaches as well
# point to point nodes for all shapefiles might be the simplest

def get_indexes(row, col, radius):
    raise NotImplementedError

def get_tranmissivity(row, col, layer, radius=5000):
    """
    get the transmissivity as an average of teh area
    :param row: well row
    :param col: well column
    :param layer: well layer
    :param radius: radius to consider (m) this is considered the square rounded number of cells that the radius
                   encounters
    :return:
    """
    #todo make a _getindexs
    raise NotImplementedError

def get_aquitard_kv(row, column, layer, radius=5000):
    """
    get the transmissivity as an average of teh area
    :param row: well row
    :param col: well column
    :param layer: well layer
    :param radius: radius to consider (m) this is considered the square rounded number of cells that the radius
                   encounters
    :return:
    """
    raise NotImplementedError

def get_aquitard_sat_thickness(row, column, layer, radius=5000):
    """
    get the transmissivity as an average of teh area
    :param row: well row
    :param col: well column
    :param layer: well layer
    :param radius: radius to consider (m) this is considered the square rounded number of cells that the radius
                   encounters
    :return:
    """
    raise NotImplementedError


def get_separation_distance():
    raise NotImplementedError