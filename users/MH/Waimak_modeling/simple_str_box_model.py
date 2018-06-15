"""
Author: matth
Date Created: 12/06/2017 2:59 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from supporting_data_path import sdp
import flopy_mh as flopy
import os

model_ws = r"C:\Users\MattH\Downloads\test_box_model"
if not os.path.exists(model_ws):
    os.makedirs(model_ws)

m = flopy.modflow.Modflow('test',"{}/models_exes/mf2005.exe".format(sdp),model_ws=model_ws)
bas = flopy.modflow.ModflowBas(m,0)
lpf = flopy.modflow.ModflowLpf(m,)