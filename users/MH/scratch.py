# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 8:37 AM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.percentage_reduction_maps import get_pa_reductions


if __name__ == '__main__':
    t = get_pa_reductions()
    print(t)
