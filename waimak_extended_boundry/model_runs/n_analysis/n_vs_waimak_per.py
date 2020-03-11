# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/03/2018 3:16 PM
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
from other_functions.stats import LR
from matplotlib.offsetbox import AnchoredText
import os
import numpy as np
from waimak_extended_boundry.model_runs.n_analysis.nitrate_at_key_receptors import get_well_ids


def n_vs_wai_streams(outdir):
    # saltwater increases with increasing alpine river water because the Northern boundary flux was mistakenly not set
    # as LSR in the EMMA models.
    npath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\raw_stocastic_set_str_data.csv"
    waipath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_str_data.csv"

    if not os.path.exists(os.path.join(outdir, 'plots')):
        os.makedirs(os.path.join(outdir, 'plots'))

    n = pd.read_csv(npath, index_col=0)
    sites = n.keys()
    wai = pd.read_csv(waipath, index_col=0)
    data = pd.merge(n, wai, left_index=True, right_index=True, suffixes=('_n', '_wai'))
    outdata = pd.DataFrame(index=sites, columns=['m', 'b', 'adj_r2'])
    for site in sites:
        if site == 'cust_skewbridge':
            n_temp = data.loc[:, '{}_n'.format(site)]
            all_aloss = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\ashley_losses_and_chch_con.csv", index_col=0)
            # ashley losses in cumics (pulled out from smc_visualisation\upper_ashley_loss)
            wai_temp = all_aloss.loc[n_temp.index,'ash']
        elif data.loc[:, '{}_wai'.format(site)].isnull().all():
            continue
        else:
            n_temp = data.loc[:, '{}_n'.format(site)]
            wai_temp = data.loc[:, '{}_wai'.format(site)]
        fig, ax = plt.subplots(figsize=(18.5, 9.5))
        model = LR(wai_temp, n_temp)
        ax.scatter(wai_temp, n_temp)
        ax.plot(wai_temp, model.predict(wai_temp))
        anchored_text = AnchoredText("formula: {}\nadj_R2: {}".format(model.formula, model.adj_rval), loc=2)
        ax.add_artist(anchored_text)
        ax.set_ylabel('N')
        ax.set_xlabel('alpine_river_fraction (or ashley_loss)')
        ax.set_title(site)
        fig.savefig(os.path.join(outdir, 'plots', '{}.png'.format(site)))
        plt.close(fig)
        outdata.loc[site, 'm'] = model.slope
        outdata.loc[site, 'b'] = model.intercept
        outdata.loc[site, 'adj_r2'] = model.adj_rval
    outdata.to_csv(os.path.join(outdir,'stream_regression_data.csv'))

def n_vs_wai_wells(outdir):
    npath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\raw_stocastic_set_well_data.csv"
    waipath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_well_data.csv"

    if not os.path.exists(os.path.join(outdir, 'plots')):
        os.makedirs(os.path.join(outdir, 'plots'))

    n = pd.read_csv(npath, index_col=0)
    wai = pd.read_csv(waipath, index_col=0)
    data = pd.merge(n, wai, left_index=True, right_index=True, suffixes=('_n', '_wai'))

    # sort out sites
    wells = get_well_ids()
    wells.loc[wells.Zone.notnull(),'Zone'] = ['wdc_{}'.format(e) for e in wells.loc[wells.Zone.notnull(),'Zone']]
    zone_sets = [set(wells.Zone[wells.Zone.notnull()]),
            set(wells.Zone_1[wells.Zone_1.notnull()]),
            set(wells.zone_2[wells.zone_2.notnull()])]
    sites = list(zone_sets[0] | zone_sets[1] | zone_sets[2])

    outdata = pd.DataFrame(index=sites, columns=['m', 'b', 'adj_r2'])
    for zone_set, key in zip(zone_sets, ['Zone','Zone_1','zone_2']):
        for site in zone_set:
            idxs = wells.loc[wells[key]==site].index

            n_temp = data.loc[:, ['{}_n'.format(e)for e in idxs]]
            n_temp = np.nanmean(n_temp.values, axis=1)
            wai_temp = data.loc[:, ['{}_wai'.format(e)for e in idxs]]
            wai_temp = np.nanmean(wai_temp.values, axis=1)

            fig, ax = plt.subplots(figsize=(18.5, 9.5))
            model = LR(wai_temp, n_temp)
            ax.scatter(wai_temp, n_temp)
            ax.plot(wai_temp, model.predict(wai_temp))
            anchored_text = AnchoredText("formula: {}\nadj_R2: {}".format(model.formula, model.adj_rval), loc=2)
            ax.add_artist(anchored_text)
            ax.set_ylabel('N')
            ax.set_xlabel('alpine_river_fraction')
            ax.set_title(site)
            fig.savefig(os.path.join(outdir, 'plots', '{}.png'.format(site)))
            plt.close(fig)
            outdata.loc[site, 'm'] = model.slope
            outdata.loc[site, 'b'] = model.intercept
            outdata.loc[site, 'adj_r2'] = model.adj_rval
    outdata.to_csv(os.path.join(outdir,'well_regression_data.csv'))

if __name__ == '__main__':
    base_outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_vs_wai_regressions"
    n_vs_wai_streams(os.path.join(base_outdir,'streams'))
    n_vs_wai_wells(os.path.join(base_outdir,'wells'))
    print('done')

