# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 8:37 AM
"""
from __future__ import division
from core import env
import os
import tempfile
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.percentage_reduction_maps import get_pa_reductions
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import flopy_mh as flopy


if __name__ == '__main__':
    m = flopy.modflow.Modflow.load(r"C:\Users\MattH\Desktop\can probably delete\AshOpt_simple_modpath\AshOpt_modpath_base.nam")
    print('done')