# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/05/2018 8:47 AM
"""

from __future__ import division
from core import env
from setup_run_models import setup_run_mt3d_suite, extract_data, pc5_ftl_repo
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.mt3d_wrapper import \
    get_sft_stress_period_data, get_ssm_stress_period_data, setup_run_mt3d
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import os
from warnings import warn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from core.stats.LR_class import LR
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.interzone_n import get_interzone_n, \
    _np_describe
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.nitrate_at_key_receptors import \
    get_n_at_points_nc, get_str_ids, get_well_ids

alpine_fractions = {
    # site name: modified n value

    # WDC wells
    # ## leave as is:
    # ###### fernside
    # ###### manderville
    # ###### oxford urban
    # ###### pegasus
    # ###### poyntzs rd
    # ###### waikuku

    # cust emmma calculation s median of 18% alpine
    'wdc_Cust': .18,

    # kaiapoi # from emma medain of 32% alpine
    'wdc_Kaiapoi': .32,

    # kairaki  #from emma medain from kaiakpoi of 32% alpine
    'wdc_Kairaki': .32,

    # ohoka # from EMMA median of 20%
    'wdc_Ohoka': .20,

    # rangiora #from emma medain from kaiakpoi of 32% alpine
    'wdc_Rangiora': .32,

    # west eyreton as west eyreton wells  # from emma 19% alpine
    'wdc_West Eyreton': .19,

    # private wells
    # ## leave as is:
    # 'Fernside',
    # 'Flaxton',
    # 'Horellville',
    # 'Mandeville',
    # 'North East Eyrewell_shallow',
    # 'North West Eyrewell_shallow',
    # 'Rangiora',
    # 'Swannanoa_shallow',
    # 'Waikuku',
    # 'Woodend - Tuahiwi',
    # 'West Eyreton_shallow',


    # Clarkville , scale similar to kaiapoi at island road waimakariri component 20% alpine river water
    'Clarkville': .20,

    # Cust, similar to cust WDC well emmma calculation s median of 18% alpine
    'Cust': .18,

    # 'Eyreton_deep', use our EMMA data 32% alpine river water
    'Eyreton_deep': .32,

    # 'Eyreton_shallow', as kaiapoi at island road component 20% alpine river water
    'Eyreton_shallow': .20,

    # 'North East Eyrewell_deep', # look for EMMA lump with NW eyrewell 28% alpine river water from emma
    'North East Eyrewell_deep': .28,

    # 'North West Eyrewell_deep', # look for EMMA lump with NE eyrewell 28% alpine river water from emma
    'North West Eyrewell_deep': .28,

    # 'Ohoka_deep', # as ohoka WDC well  from EMMA median of 20%
    'Ohoka_deep': .20,

    # 'Ohoka_shallow', # as ohoka stream # from emma we can expect a median of 12% waimak water
    'Ohoka_shallow': .12,

    # 'Springbank', #as cust (and cust wdc) emmma calculation s median of 18% alpine
    'Springbank': .18,

    # 'Summerhill', # look for emma 26% alpine river water
    'Summerhill': .26,

    # 'Swannanoa_deep', as WDC west eyreton # from emma 19% alpine
    'Swannanoa_deep': .19,

    # 'West Eyreton_deep', as swannanowa/WDC west eyreton # from emma 19% alpine
    'West Eyreton_deep': .19,
    # site name: modified n value

    # not modifing
    # cam_bramley_s: no change, we do not have EMMA data and it is unlikely to lead to any real changes (e.g. NOF bands)
    # cam at mashes as cam at bramleys
    # cust_skewbridge: no change, it is unlikely to lead to any real changes (e.g. NOF bands)
    # northbrook and southbrook at marshes: no change, it is unlikely to lead to any real changes (e.g. NOF bands)

    # the below modified from regressions created from n_vs_waimak_per values are in mg/l
    # courtaney at the kaiapoi
    'courtenay_kaiapoi_s': 0.8,  # from emma we expect a median of 8% waimak water

    # kaiapoi (silver stream) at harpers road
    'kaiapoi_harpers_s': .23,  # from emma we can expect a median of 23% waimak water

    # kaiapoi at island road
    'kaiapoi_island_s': .20,  # from emma we can expect a median of 20% waimak water

    # ohoka at island road
    'ohoka_island_s': .12  # from emma we can expect a median of 12% waimak water
}
sites = {
    # streams
    'cust_skewbridge',
    'cam_bramleys_s',
    'cam_marshes_s',
    'courtenay_kaiapoi_s',
    'kaiapoi_harpers_s',
    'kaiapoi_island_s',
    'northbrook_marshes_s',
    'ohoka_island_s',
    'southbrook_marshes_s',

    # wdc_wells
    'wdc_Kairaki',
    'wdc_Oxford Urban',
    'wdc_Ohoka',
    'wdc_Fernside',
    'wdc_Rangiora',
    'wdc_Mandeville',
    'wdc_West Eyreton',
    'wdc_Kaiapoi',
    'wdc_Waikuku',
    'wdc_Poyntzs Road',
    'wdc_Pegasus',
    'wdc_Cust',

    # private wells
    'Fernside',
    'Flaxton',
    'Horellville',
    'Mandeville',
    'North East Eyrewell_shallow',
    'North West Eyrewell_shallow',
    'Rangiora',
    'Swannanoa_shallow',
    'Waikuku',
    'Woodend - Tuahiwi',
    'West Eyreton_shallow',
    'Clarkville',
    'Cust',
    'Eyreton_deep',
    'Eyreton_shallow',
    'North East Eyrewell_deep',
    'North West Eyrewell_deep',
    'Ohoka_deep',
    'Ohoka_shallow',
    'Springbank',
    'Summerhill',
    'Swannanoa_deep',
    'West Eyreton_deep',

}


def get_well_nums_for_group():
    str_ids = get_str_ids()
    wells = get_well_ids()
    zone_sets = [set(wells.Zone[wells.Zone.notnull()]),
                 set(wells.Zone_1[wells.Zone_1.notnull()]),
                 set(wells.zone_2[wells.zone_2.notnull()])]
    temp = {}
    for zone_set, key in zip(zone_sets, ['Zone', 'Zone_1', 'zone_2']):
        for zone in zone_set:
            if key == 'Zone':
                use_zone = 'wdc_{}'.format(zone)
            else:
                use_zone = zone
            temp[use_zone] = wells.loc[wells[key] == zone].index.values

    outdata = {}
    for key in sites:
        if key in str_ids:
            continue
        outdata[key] = temp[key]

    return outdata


def setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, ftl_repo=pc5_ftl_repo, dt0=1e4, ttsmax=1e5):
    ssm = get_ssm_stress_period_data()
    sft = get_sft_stress_period_data()

    # get scon from a single run to try and speed up mt3d runs
    scon = None
    # in future I could probably speed up runs by running a single mt3d run.  a test showed a run time of
    # 1412s for a scon of zero and 1335 s for a scon of another run... doesn't seem worthwhile


    setup_run_mt3d_suite(base_mt3d_dir, ftl_repo,
                         ssm_crch={0: rch_con},
                         ssm_stress_period_data={0: ssm},
                         sft_spd={0: sft}, dt0=dt0, ttsmax=ttsmax, scon=scon)

    extract_data(base_mt3d_dir, outfile=out_nc, description=nc_description, nname='mednload', units='g/m3')


def extract_receptor_data(scenario_paths, cbc_paths, outdir):  # todo check/test
    """

    :param scenario_paths: dictionary scenario:path to nc file;  scenario must have gmp or cmp in the scenario name
    :param cbc_paths: either the path to the cbc_netcdf (string) or a dictionary scenario and path to teh cbc_netcdf
    :return:
    """
    # extract the raw data from the model runs for all receptors intra and interzone...

    # interzone receptors including n load uncertainty
    interzone_outdir = os.path.join(outdir, 'interzone')
    if not os.path.exists(interzone_outdir):
        os.makedirs(interzone_outdir)
    get_interzone_n(scenario_paths.keys(),
                    os.path.join(interzone_outdir, 'all_n_interzone.csv'))
    data = pd.read_csv(os.path.join(interzone_outdir, 'all_n_interzone.csv'),
                       index_col=[0, 1], header=[0, 1])
    data = data.reorder_levels(['stat', 'scenario'], axis=1)
    data['1%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_1ths.csv'))
    data['5%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_5ths.csv'))
    data['10%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_10ths.csv'))
    data['25%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_25ths.csv'))
    data['50%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_50ths.csv'))
    data['75%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_75ths.csv'))
    data['95%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_95ths.csv'))
    data['99%'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_99ths.csv'))
    data['mean'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_mean.csv'))
    data['std'].to_csv(os.path.join(interzone_outdir, 'n_data_interzone_std.csv'))

    # intrazone data
    intrazone_dir = os.path.join(outdir, 'waimakariri_zone')

    # raw data
    raw_dir = os.path.join(intrazone_dir, 'raw_model_data')
    corrected_dir = os.path.join(intrazone_dir, 'corrected_model_data')
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)
    if not os.path.exists(corrected_dir):
        os.makedirs(corrected_dir)
    scenarios = scenario_paths.keys()
    pers_names = _np_describe(smt.get_empty_model_grid()).index.values
    outdata = pd.DataFrame(index=sites,
                           columns=pd.MultiIndex.from_product((scenarios, pers_names),
                                                              names=['scenario', 'stat']), dtype=float)
    well_sites = get_well_nums_for_group()
    cmp_waimak_data = {
        'stream': pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations an"
                                      r"d results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set"
                                      r"_str_data.csv"), index_col=0).transpose(),
        'well': pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations a"
                                    r"nd results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_"
                                    r"well_data.csv"),index_col=0).transpose()
    }
    for scen, path in scenario_paths.items():
        plot_dir = os.path.join(corrected_dir,'{}_plots'.format(scen))
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        if isinstance(cbc_paths, str):
            use_cbc_path = cbc_paths
        elif isinstance(cbc_paths, dict):
            use_cbc_path = cbc_paths[scen]
        else:
            raise ValueError('unexpected types for cbc_paths: {}'.format(type(cbc_paths)))
        get_n_at_points_nc(outdir, nsmc_nums='all', ucn_var_name='mednload',
                           ucn_nc_path=path,
                           cbc_nc_path=use_cbc_path,
                           missing_str_obs='raise')

        all_well_data = pd.read_csv(os.path.join(raw_dir, 'raw_stocastic_set_well_data.csv'), index_col=0).transpose()
        all_str_data = pd.read_csv(os.path.join(raw_dir, 'raw_stocastic_set_str_data.csv'), index_col=0).transpose()
        for site in sites:
            if site in get_str_ids():
                recp_type = 'stream'
                # do stream stuff
                raw_data = all_str_data.loc[site]
            else:
                # do well stuff
                recp_type = 'well'
                raw_data = all_well_data.loc[well_sites[site]]

            if site in alpine_fractions.keys():
                # correct the data for EMMA concerns
                raw_data = correct_alpine_river(site, cmp_waimak_data, raw_data, well_sites, plot_dir)

            # add the stocastic N load
            stocastic_data = add_stocastic_load(site, raw_data)
            for per in pers_names:
                outdata.loc[site, (scen, per)] = stocastic_data[per]

    outdata.to_csv(os.path.join(corrected_dir, 'all_n_waimak_zone.csv'))
    outdata = outdata.reorder_levels(['stat', 'scenario'], axis=1)
    outdata['1%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_1ths.csv'))
    outdata['5%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_5ths.csv'))
    outdata['10%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_10ths.csv'))
    outdata['25%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_25ths.csv'))
    outdata['50%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_50ths.csv'))
    outdata['75%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_75ths.csv'))
    outdata['95%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_95ths.csv'))
    outdata['99%'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_99ths.csv'))
    outdata['mean'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_mean.csv'))
    outdata['std'].to_csv(os.path.join(interzone_outdir, 'n_data_waimak_zone_std.csv'))


def add_stocastic_load(site, base_data):
    base_data = np.atleast_1d(base_data)
    if site in get_str_ids():
        # stream
        base_dir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and "
                           r"results\ex_bd_va\n_results\nwaimak_springfeds\stocastic_n_nwaimak_springfeds\without_"
                           r"trans\nconc_cmp_second_tranche\raw_data")  # cmp as gmp has a reduction associated
    elif 'wdc' in site:
        # wdc_well
        base_dir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                           r"s\ex_bd_va\n_results\wdc_use_mix\stocastic_n_wdc_use_mix\without_trans\nconc_cmp_wdc_use"
                           r"_mix\raw_data")  # cmp as gmp has a reduction associated
    else:
        # private well
        base_dir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_"
                           r"bd_va\n_results\private_wells_90\stocastic_n_private_wells_90\without_trans\nconc_cmp_"
                           r"private_wells_90_named_right\raw_data")  # cmp as gmp has a reduction associated
    modifiers = np.loadtxt(os.path.join(base_dir, '{}.txt'.format(site.replace('wdc_', ''))))

    raw_data = base_data.copy().flatten()
    all_n = raw_data[:, np.newaxis] * modifiers[np.newaxis, :]
    outdata = _np_describe(all_n)
    return outdata

def correct_alpine_river(site, waimak_data, n_data,well_sites, plot_dir):
    warn('correcting data for EMMA results please review plots in: {}'.format(plot_dir))
    if site in get_str_ids():
        # stream stuff
        wai = waimak_data['stream'].loc[site]
        n = n_data.loc[site]
    else:
        wai = waimak_data['well'].loc[well_sites[site]].mean(axis=0)
        n = n_data.loc[well_sites[site]].mean(axis=0)
        # well stuff
    plot_data = pd.merge(pd.DataFrame(wai,columns=['wai']),pd.DataFrame(n,columns=['n']),
                         left_index=True, right_index=True)
    model = LR(plot_data.wai, plot_data.n)
    outdata = model.predict(alpine_fractions[site])
    n_temp = plot_data.n
    wai_temp = plot_data.wai
    fig, ax = plt.subplots(figsize=(18.5, 9.5))
    model = LR(wai_temp, n_temp)
    ax.scatter(wai_temp, n_temp)
    ax.plot(wai_temp, model.predict(wai_temp))
    anchored_text = AnchoredText("formula: {}\nadj_R2: {}".format(model.formula, model.adj_rval), loc=2)
    ax.add_artist(anchored_text)
    ax.set_ylabel('N')
    ax.set_xlabel('alpine_river_fraction (or ashley_loss)')
    ax.set_title(site)
    fig.savefig(os.path.join(plot_dir, 'plots', '{}.png'.format(site)))
    plt.close(fig)

    return outdata

if __name__ == '__main__':
    print('done')