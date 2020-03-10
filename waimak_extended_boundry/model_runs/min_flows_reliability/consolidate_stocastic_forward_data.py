# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/03/2018 1:55 PM
"""

from __future__ import division
import env
import pandas as pd
import numpy as np
import os
from glob import glob
from waimak_extended_boundry.model_run_tools import \
    get_samp_points_df


# NsmcReal002160_non_cc_forward_runs_2018-03-24_results\NsmcReal002160_relative_data.csv"

def consolidate_forward_runs(base_dir, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    streams = get_samp_points_df()
    streams = streams.loc[streams.m_type == 'min_flow'].index

    # get nsmc_nums where the current converged (e.g. no exception raised in process)
    rel_data_paths = glob(os.path.join(base_dir, "*_results\*_relative_data.csv"))
    meta_data_paths = [e.replace('_relative_data.csv', '_meta_data.csv') for e in rel_data_paths]
    nsmc_nums = [int(os.path.basename(e).split('_')[0].replace('NsmcReal', '')) for e in rel_data_paths]

    scenarios = list(pd.read_csv(meta_data_paths[0]).iloc[:, 0])

    # initalise outdatas
    outdatas = {}
    for site in streams:
        outdatas[site] = pd.DataFrame(index=list(set(nsmc_nums)), columns=scenarios)

    # add only the data from those that converged
    for nsmc_num, rel_path, meta_path in zip(nsmc_nums, rel_data_paths, meta_data_paths):
        meta_data = pd.read_csv(meta_path, index_col=0)
        rel_data = pd.read_csv(rel_path, skiprows=1, index_col=0)
        scens_converged = list(set(meta_data.index[meta_data.converged]) - {'current'})
        for site in streams:
            outdatas[site].loc[nsmc_num, scens_converged] = rel_data.loc[site, scens_converged]

    pers = [0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95, 0.99]
    if not os.path.exists(os.path.join(outdir, 'raw_data')):
        os.makedirs(os.path.join(outdir, 'raw_data'))

    for site in streams:
        outdatas[site] = outdatas[site].astype(float)
        outdatas[site][outdatas[site]<-10] = np.nan
        outdatas[site].to_csv(os.path.join(outdir, 'raw_data', '{}_raw_data.csv'.format(site)))
        outdatas[site].astype(float).describe(percentiles=pers).transpose().to_csv(
            os.path.join(outdir, '{}_summary.csv'.format(site)))
    consolidate_by_scenario(outdir)


def consolidate_by_scenario(base_dir):
    scenarios = ['current_w_ncar',
                 'full_abs',
                 'full_abs_allo',
                 'full_allo_cur_use',
                 'mod_period',
                 'mod_period_w_ncar',
                 'naturalised',
                 'pc5_80',
                 'pc5_80_full_allo_cur_usage',
                 'pc5_80_wil_eff',
                 'pc5_no_pump_reduc',
                 'pc5_no_pump_reduc_wil_eff',
                 'super_gmp',
                 'wil_eff']

    metrics = ['count',
               'mean',
               'std',
               'min',
               '1%',
               '5%',
               '10%',
               '25%',
               '50%',
               '75%',
               '90%',
               '95%',
               '99%',
               'max']
    streams = get_samp_points_df()
    streams = streams.loc[streams.m_type == 'min_flow'].index
    outdatas = {}
    for scen in scenarios:
        outdatas[scen] = pd.DataFrame(index=streams, columns=metrics)

    for site in streams:
        data = pd.read_csv(os.path.join(base_dir, '{}_summary.csv'.format(site)), index_col=0)
        for scen in scenarios:
            outdatas[scen].loc[site,metrics] = data.loc[scen,metrics]

    for scen in scenarios:
        outdatas[scen].to_csv(os.path.join(base_dir, '{}_summary.csv'.format(scen)))




if __name__ == '__main__':
    consolidate_forward_runs(env.gw_met_data(r"mh_modeling\stocastic_forward"),
                             env.gw_met_data(r"mh_modeling\stocastic_forward\condensed_data"))
