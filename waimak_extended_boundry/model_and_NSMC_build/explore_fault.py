# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 22/07/2017 9:12 PM
"""

from __future__ import division
from waimak_extended_boundry.model_and_NSMC_build.m_packages.sfr2_packages import _get_reach_data
from waimak_extended_boundry.model_and_NSMC_build.m_packages.wel_packages import get_wel_spd
from waimak_extended_boundry.model_and_NSMC_build.m_packages.drn_packages import _get_drn_spd
from waimak_extended_boundry import smt
import pandas as pd
import numpy as np


def get_all_bcs():
    sfr_data = pd.DataFrame(_get_reach_data(smt.reach_v))
    sfr_data['bc_type'] = 'sfr'

    well_data = get_wel_spd(smt.wel_version).rename(columns={'row':'i','col':'j', 'layer':'k'})
    well_data['bc_type'] = 'well'

    drn_data = _get_drn_spd(smt.reach_v,smt.wel_version)
    drn_data['bc_type'] = 'drn'
    drn_data = drn_data.loc[~drn_data.group.str.contains('carpet')] # don't worry about carpet drains

    all_data = pd.concat((sfr_data,well_data,drn_data)).reset_index()

    for i in all_data.index:
        row, col = all_data.loc[i,['i','j']]
        x,y = smt.convert_matrix_to_coords(row,col)
        all_data.loc[i,'mx'] = x
        all_data.loc[i,'my'] = y
    return all_data

def get_fault_zone():
    fault_zone = np.repeat(smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/Pegusus_bay_fault.shp".format(smt.sdp),
                                                         'Id',True)[np.newaxis,:,:],smt.layers,axis=0)
    no_flow = smt.get_no_flow()
    fault_zone[np.isclose(no_flow,0)] = np.nan

    fault_zone+=1
    fault_zone[np.isnan(fault_zone)] = 0
    return fault_zone

if __name__ == '__main__':
    test = get_fault_zone()
    print('done')