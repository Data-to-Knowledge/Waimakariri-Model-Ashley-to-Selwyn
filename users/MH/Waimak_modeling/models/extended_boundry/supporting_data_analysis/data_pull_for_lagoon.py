# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/06/2018 2:16 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import get_model
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import pandas as pd

if __name__ == '__main__':
    m = get_model('NsmcBase')
    #row,col = smt.convert_coords_to_matix(1575962,5203005) # for lagoon
    row,col = smt.convert_coords_to_matix(1556717,5193196) # for infiltration test #todo
    elv_db = smt.calc_elv_db()

    outdata = pd.DataFrame(index=range(11))
    outdata.loc[:,'kh'] = m.upw.hk.array[:,row,col]
    outdata.loc[:,'kv'] = m.upw.vka.array[:,row,col]
    outdata.loc[:,'thickness'] = (elv_db[0:-1] - elv_db[1:])[:,row,col]
    #outdata.to_csv(r"C:\Users\MattH\Downloads\data_for_lagoon.csv") # for lagoon
    outdata.to_csv(r"T:\Temp\temp_gw_files\data_for_itest.csv") # for infiltration test