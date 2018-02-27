# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/08/2017 2:17 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np


def gain_loss_flows_target_uncert(data):
    """
    calculate the gain and loss target for the data
    :param data: n1,n2,R1,R2,Y1,Y2 for each site 1 is upstream 2 is downstream for
    :return:
    """
    n1, n2, R1, R2, y1, y2 = data.loc[:['n1','n2','R1','R2','y1','y2']]
    y1_error = 0.085*y1 + 1.96*(n1 / (n1 - 2) * y1**2 * (1-R1)*(1+1/n1))**0.5
    y2_error = 0.085*y2 + 1.96*(n2 / (n2 - 2) * y2**2 * (1-R2)*(1+1/n2))**0.5

    y1_error.loc[pd.isnull(y1_error)] = 0
    y2_error.loc[pd.isnull(y2_error)] = 0


    data.loc[:,'min_uncert'] = (y1_error**2 + y2_error**2)**0.5

    return data

