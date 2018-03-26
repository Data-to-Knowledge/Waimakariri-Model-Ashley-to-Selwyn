# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/03/2018 1:55 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
import os
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction import \
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
        outdatas[site] = pd.DataFrame(index=nsmc_nums, columns=scenarios)

    # add only the data from those that converged
    for nsmc_num, rel_path, meta_path in zip(nsmc_nums, rel_data_paths, meta_data_paths):
        meta_data = pd.read_csv(meta_path, index_col=0)
        rel_data = pd.read_csv(rel_path, skiprows=1, index_col=0)
        scens_converged = list(set(meta_data.index[meta_data.converged])-{'current'})
        for site in streams:
            outdatas[site].loc[nsmc_num, scens_converged] = rel_data.loc[site, scens_converged]

    pers = [0.01,0.05,0.1,0.25,0.50,0.75,0.9,0.95,0.99]
    for site in streams:
        outdatas[site].to_csv(os.path.join(outdir,'raw_data', '{}_raw_data.csv'.format(site)))
        outdatas[site].describe(percentiles=pers).transpose().to_csv(os.path.join(outdir, '{}_summary.csv'.format(site)))


if __name__ == '__main__':
    consolidate_forward_runs(env.gw_met_data(r"mh_modeling\stocastic_forward"),
                             env.gw_met_data(r"mh_modeling\stocastic_forward\condensed_data"))
