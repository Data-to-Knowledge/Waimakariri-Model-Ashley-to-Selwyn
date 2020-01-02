"""
Author: matth
Date Created: 3/04/2017 3:38 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np

def read_visas_matrix (path):
    """
    read the vistas matrix to a np array size 190 x 365L
    :param path:
    :return:
    """
    temp = pd.read_table(path, delim_whitespace=True,header=None)
    outdata = np.array(temp)

    return outdata

def write_vistas_matrix (path, array):
    np.savetxt(path, array)

bc_path = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/supporting_data_for_scripts/BCs.txt")

no_flow = np.invert((np.loadtxt(bc_path)[0:190,:]).astype(bool)) # needs updating for new model

if __name__ == '__main__':
    test = read_visas_matrix(r"C:\Users\MattH\Downloads\BHMODEL_OPTFLOW2.DAT")
    #print(no_flow.shape)