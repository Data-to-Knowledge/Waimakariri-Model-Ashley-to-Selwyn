# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/10/2017 7:40 PM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import get_wel_spd
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.drn_packages import _get_drn_spd
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.sfr2_packages import _get_reach_data
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

if __name__ == '__main__':
    well_data = get_wel_spd(3)
    well_data['bc_type'] = 'well'
    well_data = well_data.loc[well_data.layer==0]
    well_data = well_data.rename(columns={'row':'i','col':'j','layer':'k'})
    drn_data = _get_drn_spd(1,1)
    carpet_names = [
        'ash_carpet',
        'chch_carpet',
        'cust_carpet',
        'down_lincoln',
        'down_selwyn',
        'up_lincoln',
        'up_selwyn',
    ]
    drn_data = drn_data.loc[~np.in1d(drn_data.group,carpet_names)]
    drn_data['bc_type'] = 'drn'
    drn_data = smt.add_mxmy_to_df(drn_data)
    sfr_data = pd.DataFrame(_get_reach_data(1))
    sfr_data['bc_type'] = 'sfr'
    sfr_data = smt.add_mxmy_to_df(sfr_data)

    all_data = pd.concat((well_data,drn_data,sfr_data))

    all_data['dup'] = all_data.duplicated(subset=['i','j','k'], keep=False)
    all_data.to_csv(r"C:\Users\MattH\Downloads\check_overlap.csv")
