"""
Author: matth
Date Created: 12/05/2017 2:13 PM
"""

from __future__ import division
from core import env
import flopy
import numpy as np
import pandas as pd
from well_values import get_all_well_data,get_well_postions
from m_wraps.base_modflow_wrapper import get_base_mf_ss
from drain_concentration import get_drn_concentration, get_drn_samp_pts_dict, get_drn_samp_pts
from users.MH.Waimak_modeling.supporting_data_path import sdp

def concentration_at_wells(well_nums, con, mt3d_kskper = None):
    well_nums = np.atleast_1d(well_nums)
    outdata = pd.Series(index=well_nums)
    well_data = get_all_well_data()
    if isinstance(con, str):
        con = flopy.utils.UcnFile(con)
    if mt3d_kskper is None:
        mt3d_kskper = con.get_kstpkper()[-1] # use the last time step is none
    gw_conc_data = con.get_data(mt3d_kskper)

    for well in well_nums:
        if well not in list(well_data.index):
            layer, row, col = get_well_postions(well)
        else:
            layer = well_data.loc[well, 'layer']
            row = well_data.loc[well, 'row']
            col = well_data.loc[well, 'col']
        outdata.loc[well] = gw_conc_data[layer, row, col]

    return outdata

def get_n_at_targets(con,m=None, mt3d_kskper=None, mf_kskper=None): #add multipel time periods or all some day?
    if isinstance(con, str):
        con = flopy.utils.UcnFile(con)

    #wells
    # load target wells
    target_wells = pd.read_csv("{}/inputs/wells/N calibration targets for GW model_from_fouad.csv".format(sdp))
    target_wells = target_wells.set_index('WELL_NO')
    target_wells['target_g_3'] = target_wells['N_kg_m3']*1000

    target_wells['model_g_m3'] = concentration_at_wells(target_wells.index,con,mt3d_kskper)

    #drains
    if m is None:
        m = get_base_mf_ss()
    drn_points = get_drn_samp_pts()
    drain_values = get_drn_concentration(drn_points, m, con, mt3d_kskper, mf_kskper)

    #streams #not done

    return target_wells, drain_values


if __name__ == '__main__':
    wells, drains = get_n_at_targets(r"D:\mattH\python_wm_runs\first_mt3d2\cmp_layer\MT3D001.UCN")
    wells.to_csv(r"C:\Users\matth\Desktop\well_targets.csv")
    drains.to_csv(r"C:\Users\matth\Desktop\drain_targets.csv")


