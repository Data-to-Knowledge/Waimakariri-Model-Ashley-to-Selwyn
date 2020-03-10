"""
Author: matth
Date Created: 2/05/2017 4:22 PM
"""

from __future__ import division

from basic_tools import calc_elv_db
from drain_concentration import get_drn_samp_pts_dict
from m_wraps.base_modflow_wrapper import get_base_mf_ss
from users.MH.Waimak_modeling.model_runs.inital_runs.mpath_runs.flow_paths_wel import get_domestic_particles
from well_values import get_race_data, get_original_well_data, get_all_well_data, get_model_well_full_consented
from users.MH.Waimak_modeling.model_runs.stream_depletion_modeling.transient_inputs import get_mean_month_well_scaling, \
    get_monthly_str_scale
from get_str_rch_values import get_base_str, get_base_rch
from users.MH.Waimak_modeling.model_runs.stream_depletion_modeling.s_cal_targets import get_s_cal_data

def recalc_all_pick_base():
    get_race_data(True)
    get_original_well_data(True)
    get_drn_samp_pts_dict(True)
    calc_elv_db(True)
    get_all_well_data(True)
    get_base_mf_ss(recalc=True)
    get_base_str(True)
    get_base_rch(True)
    get_model_well_full_consented(True)

def recalc_all_pick_runs():
    get_domestic_particles(True)
    get_mean_month_well_scaling(True)
    get_monthly_str_scale(True)
    get_s_cal_data(True)


if __name__ == '__main__':
    run_type = 3
    if run_type == 1:
        recalc_all_pick_base()
    elif run_type == 2:
        recalc_all_pick_runs()
    elif run_type == 3:
        recalc_all_pick_base()
        recalc_all_pick_runs()
    elif run_type == 4:
        get_base_mf_ss(recalc=True)