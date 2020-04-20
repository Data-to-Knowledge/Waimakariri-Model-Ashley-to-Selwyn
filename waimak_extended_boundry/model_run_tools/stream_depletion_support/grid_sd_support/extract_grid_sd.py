# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/10/2017 10:20 AM
"""

from __future__ import division
import flopy_mh as flopy
import pandas as pd
import glob
import numpy as np
import os
from waimak_extended_boundry.model_run_tools.metadata_managment.convergance_check import modflow_converged
from waimak_extended_boundry.model_run_tools.data_extraction.data_from_streams import get_flux_at_points, get_samp_points_df
from waimak_extended_boundry.model_run_tools.stream_depletion_support.grid_sd_support.ss_grid_sd_setup import grid_wells, get_base_grid_sd_path


def calc_stream_dep_grid(model_path):
    """
    calculate a dataframe for the stream depletion percentages
    :param model_path: path to the modflow model namefile with or without extension
    :return: series with index of sites defined below
    """
    model_path = model_path.replace('.nam', '')
    model_id = os.path.basename(model_path).split('_')[0]

    baseline_path = get_base_grid_sd_path(model_id)

    # id the kstpkpers and the integrater value
    kstpkpers = [(0, 0)]
    integrater = 1

    # define sites
    samp_points_df = get_samp_points_df()
    sites = list(samp_points_df[samp_points_df.m_type == 'swaz'].index)

    str_base = get_flux_at_points(sites=sites, base_path=baseline_path, kstpkpers=kstpkpers)
    str_sd = get_flux_at_points(sites=sites, base_path=model_path, kstpkpers=kstpkpers)

    outdata = (str_base.sum(axis=1) - str_sd.sum(axis=1)) * integrater * -1

    # cumulative well abstraction budget
    model_budget = flopy.utils.MfListBudget('{}.list'.format(model_path))
    temp = pd.DataFrame(model_budget.get_cumulative('WELLS_OUT'))
    abs_vol_model = temp[(temp['stress_period'] == temp['stress_period'].max()) &
                   (temp['time_step'] == temp['time_step'].max())]['WELLS_OUT'].iloc[0]

    base_budget = flopy.utils.MfListBudget('{}.list'.format(baseline_path))
    temp2 = pd.DataFrame(base_budget.get_cumulative('WELLS_OUT'))
    abs_vol_base = temp2[(temp2['stress_period'] == temp2['stress_period'].max()) &
                   (temp2['time_step'] == temp2['time_step'].max())]['WELLS_OUT'].iloc[0]
    abs_vol = abs_vol_model - abs_vol_base


    if abs_vol == 0:
        outdata *= np.nan
    elif abs_vol < 0:
        outdata *=0
        outdata += -999
    else:
        outdata *= 1 / abs_vol * 100

    if not modflow_converged('{}.list'.format(model_path)):
        outdata *= np.nan
        abs_vol *= np.nan

    return outdata, abs_vol


def calc_str_dep_all_wells_grid(out_path, base_path):
    """
    calculate the stream depletion for each well in the gridded data
    :param out_path: the path to save the csv to
    :param base_path: path to the well by well stream depletion folder
    :return: saves dataframe with index fo wells (defined from file names) and columns of sites defined in calc_stream_dep
    """
    all_paths = glob.glob('{}/*/*.nam'.format(base_path))
    all_paths = [e.replace('.nam', '') for e in all_paths]
    model_id = os.path.basename(all_paths[0]).split('_')[0]
    wells = [os.path.basename(e).replace('{}_'.format(model_id), '') for e in all_paths]
    out_path = os.path.join(os.path.dirname(out_path), '{}_{}'.format(model_id, os.path.basename(out_path)))
    outdir = os.path.dirname(out_path)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    samp_points_df = get_samp_points_df()
    sites = list(samp_points_df[samp_points_df.m_type == 'swaz'].index)

    outdata = {}
    out_per_abs_vol = {}
    for well, path in zip(wells, all_paths):
        if os.path.getsize('{}.cbc'.format(path)) == 0:  # the model fell over and nothing was written to the cbc
            per_abs_vol = np.nan
            sd = pd.Series(index=sites)
        else:
            sd, per_abs_vol = calc_stream_dep_grid(path)
        outdata[well] = sd
        out_per_abs_vol[well] = per_abs_vol
    outdata = pd.DataFrame(outdata).transpose()
    out_per_abs_vol = pd.DataFrame(out_per_abs_vol, index=['model_abs_vol']).transpose()

    outdata = pd.merge(outdata, out_per_abs_vol, right_index=True, left_index=True)

    # add additional information
    # add flux
    flux_val = float(wells[0].split('_')[-1])
    flux = grid_wells(flux_val)
    outdata = pd.merge(outdata, flux, how='left', left_index=True, right_index=True)

    # save with header
    header = (
        'SD_values for model: {}. flux: {} modflow units m and day.'
        'all sd values are relative to model_abs_vol (%); though nwt sometimes reduces a wells pumping rate if it goes dry;'
        'the flag -999 means the abstraction volume for the model was calculated at less than zero.'
        'the abstracted volume is noted in model_abs_vol(m3).  m(x;y) (e.g. modelx) and nztm(z;y) are in nztm. '
        'mid_screen_elv; mz are in m lyttleton; hor_misloc and vert_misloc are in m and vert_misloc is mz - mid_screen_elv'
        'flux is in m3/day\n'.format(model_id, flux_val))
    with open(out_path, 'w') as f:
        f.write(str(header))
    outdata.to_csv(out_path, mode='a')
