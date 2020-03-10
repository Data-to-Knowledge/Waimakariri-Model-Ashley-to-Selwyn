"""
Author: matth
Date Created: 22/05/2017 10:56 AM
"""

from __future__ import division
import numpy as np
import os
from users.MH.Waimak_modeling.supporting_data_path import sdp
import pandas as pd
from future.builtins import input

def get_ss_dist(n=10000, recalc=False):
    """
    get values from a distribution this is explicityly not set to re-run as I want values that are consistent across
    senarios
    :param n: number of values to return max 10,000
    :return: 
    """
    base_path = "{}/inputs/ss_distribution.txt".format(sdp)
    if os.path.exists(base_path) and not recalc:
        with open(base_path) as f:
          values = f.readlines()
    else:
        s_data_path = "{}/inputs/S_vs_depth Data.xlsx".format(sdp)
        sdata = pd.read_excel(s_data_path, sheetname='S for transient model')


        values = np.random.lognormal(mean=sdata.loc[1,'Mean Ln( Ss*)'], sigma=sdata.loc[1,'Sdev'],size=10000)
        values = [str(e) + '\n' for e in values]
        with open(base_path,'w') as f:
            f.writelines(values)
    values = [float(e.strip()) for e in values]

    return values[0:n]

def get_sy_dist (n=10000, recalc=False):
    """
    get values from a distribution this is explicityly not set to re-run as I want values that are consistent across
    senarios
    :param n: number of values to return max 10,000
    :return: 
    """

    base_path = "{}/inputs/sy_distribution.txt".format(sdp)
    if os.path.exists(base_path) and not recalc:
        with open(base_path) as f:
            values = f.readlines()
    else:
        s_data_path = "{}/inputs/S_vs_depth Data.xlsx".format(sdp)
        sdata = pd.read_excel(s_data_path, sheetname='S for transient model')
        values = np.random.lognormal(mean=sdata.loc[0,'Mean Ln( Ss*)'], sigma=sdata.loc[0,'Sdev'], size=10000)
        values = [str(e) + '\n' for e in values]
        with open(base_path,'w') as f:
            f.writelines(values)
    values = [float(e.strip()) for e in values]

    return values[0:n]

if __name__ == '__main__':
    cont = input("do you really want to recalculate S values y/n \n")
    if cont.lower() != 'y':
        raise ValueError('stopped to prevent recalculation')
    get_sy_dist(1,recalc=True)
    get_ss_dist(1,recalc=True)