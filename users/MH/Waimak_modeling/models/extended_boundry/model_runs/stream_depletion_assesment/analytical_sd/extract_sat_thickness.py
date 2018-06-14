# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/12/2017 8:53 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import flopy_mh as flopy
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

if __name__ == '__main__':
    elv_db = smt.calc_elv_db()
    x,y = smt.get_model_x_y()
    no_flow = smt.get_no_flow(1)
    sat_thickness = smt.get_empty_model_grid()
    hds = flopy.utils.HeadFile(r"{}\forward_supporting_models\NsmcBase_base_str_dep\NsmcBase_"
                               r"base_heads_for_str_dep.hds".format(smt.sdp)).get_data((0,0))

    idx = hds[0] > elv_db[0] #flooded cells sat thickness is the cell thickness
    sat_thickness[idx] = (elv_db[0] - elv_db[1])[idx]

    idx = hds[0] <= elv_db[0] #flooded cells sat thickness is the cell thickness
    sat_thickness[idx] = (hds[0] - elv_db[1])[idx]
    sat_thickness[sat_thickness<0] = 0

    outidx = no_flow != 0
    outdata = pd.DataFrame({'x':x[outidx],'y':y[outidx],'sat_thick': sat_thickness[outidx]})
    outdata.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_"
                   r"va\stream_depletion\saturated_thickness_layer1.csv")






