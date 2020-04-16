# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/05/2018 8:47 AM
"""

from __future__ import division
import env
from setup_run_models import setup_run_mt3d_suite, extract_data, pc5_ftl_repo
from waimak_extended_boundry.model_run_tools import get_sft_stress_period_data, get_ssm_stress_period_data, \
    get_stocastic_set
from waimak_extended_boundry import smt
import os
from warnings import warn
import numpy as np
import pandas as pd
from other_functions.stats import LR
from waimak_extended_boundry.model_run_tools.n_analysis_support.interzone_n import get_interzone_n, _np_describe
from waimak_extended_boundry.model_run_tools.n_analysis_support.nitrate_at_key_receptors import \
    get_n_at_points_nc, get_str_ids, get_well_ids

# todo it needs a touch more documenations, out of scope???

emma_dir = os.path.join(env.sdp_required, 'emma_for_n_adjustment')


def get_alpine_fractions(site=None, return_paths=False, number=1000):
    # number is the number of samples to pull out
    # if path then just send out the data (to replace alpine fraction so that I can patch this thing
    # if site is not None get teh value

    base_dir = os.path.join(emma_dir, "4_endmembers/raw_data")
    paths = {
        # site name: path to the stocastic data

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
        'wdc_Cust': ['cust-springbank_river.txt'],

        # kaiapoi # from emma medain of 32% alpine
        'wdc_Kaiapoi': [
            'kaiapoi_M35_0706_river.txt',
            'kaiapoi_M35_0788_river.txt',
            'kaiapoi_M35_0834_river.txt',
            'kaiapoi_M35_0847_river.txt',
            'kaiapoi_M35_11199_river.txt',
        ],

        # kairaki  #from emma medain from kaiakpoi of 32% alpine
        'wdc_Kairaki': [
            'kaiapoi_M35_0706_river.txt',
            'kaiapoi_M35_0788_river.txt',
            'kaiapoi_M35_0834_river.txt',
            'kaiapoi_M35_0847_river.txt',
            'kaiapoi_M35_11199_river.txt',
        ],

        # ohoka # from EMMA median of 20%
        'wdc_Ohoka': ['Ohoka-deep_river.txt'],

        # rangiora #from emma medain from kaiakpoi of 32% alpine
        'wdc_Rangiora': [
            'kaiapoi_M35_0706_river.txt',
            'kaiapoi_M35_0788_river.txt',
            'kaiapoi_M35_0834_river.txt',
            'kaiapoi_M35_0847_river.txt',
            'kaiapoi_M35_11199_river.txt',
        ],

        # west eyreton as west eyreton wells  # from emma 19% alpine
        'wdc_West Eyreton': ['swannanoa-west-eyreton_river.txt'],

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
        'Clarkville': ['kaiapoi_islands_stream_river.txt'],

        # Cust, similar to cust WDC well emmma calculation s median of 18% alpine
        'Cust': ['cust-springbank_river.txt'],

        # 'Eyreton_deep', use our EMMA data 32% alpine river water
        'Eyreton_deep': [
            'eyreton_M35_12017_river.txt',
            'eyreton_M35_12018_river.txt'
        ],

        # 'Eyreton_shallow', as kaiapoi at island road component 20% alpine river water
        'Eyreton_shallow': ['kaiapoi_islands_stream_river.txt'],

        # 'North East Eyrewell_deep', # look for EMMA lump with NW eyrewell 28% alpine river water from emma
        'North East Eyrewell_deep': ['NE-eyrewell-deep_river.txt'],

        # 'North West Eyrewell_deep', # look for EMMA lump with NE eyrewell 28% alpine river water from emma
        'North West Eyrewell_deep': ['west-NW-eyrewell_river.txt'],

        # 'Ohoka_deep', # as ohoka WDC well  from EMMA median of 20%
        'Ohoka_deep': ['Ohoka-deep_river.txt'],

        # 'Ohoka_shallow', # as ohoka stream # from emma we can expect a median of 12% waimak water
        'Ohoka_shallow': ['ohoka_islands_stream_river.txt'],

        # 'Springbank', #as cust (and cust wdc) emmma calculation s median of 18% alpine
        'Springbank': ['cust-springbank_river.txt'],

        # 'Summerhill', # look for emma 26% alpine river water
        'Summerhill': ['summerhill_river.txt'],

        # 'Swannanoa_deep', as WDC west eyreton # from emma 19% alpine
        'Swannanoa_deep': ['swannanoa-west-eyreton_river.txt'],

        # 'West Eyreton_deep', as swannanowa/WDC west eyreton # from emma 19% alpine
        'West Eyreton_deep': ['swannanoa-west-eyreton_river.txt'],
        # site name: modified n value

        # not modifing
        # cam_bramley_s: no change, we do not have EMMA data and it is unlikely to lead to any real changes (e.g. NOF bands)
        # cam at mashes as cam at bramleys
        # cust_skewbridge: no change, it is unlikely to lead to any real changes (e.g. NOF bands)
        # northbrook and southbrook at marshes: no change, it is unlikely to lead to any real changes (e.g. NOF bands)

        # the below modified from regressions created from n_vs_waimak_per values are in mg/l
        # courtaney at the kaiapoi
        # note the courtenay is possibly wrong zeb used SQ35169 which looks like a side trib
        #  when I think it should be SQ35170 which I think is the main courtenay
        # after discusstion with adrian M we decided to scrap the correction of this data as much of the saline in the
        # samples we have avalible could be derived from tidal cloride contributions
        # 'courtenay_kaiapoi_s': ['courtenay_kaiapoi_stream_river.txt'],

        # kaiapoi (silver stream) at harpers road
        'kaiapoi_harpers_s': ['kaiapoi_harpers_stream_river.txt'],
        # from emma we can expect a median of 23% waimak water

        # kaiapoi at island road
        'kaiapoi_island_s': ['kaiapoi_islands_stream_river.txt'],
        # from emma we can expect a median of 20% waimak water

        # ohoka at island road
        'ohoka_island_s': ['ohoka_islands_stream_river.txt']  # from emma we can expect a median of 12% waimak water
    }

    for key in paths.keys():
        paths[key] = [os.path.join(base_dir, e) for e in paths[key]]

    if return_paths:
        return paths

    data = []
    for path in paths[site]:
        data.append(np.loadtxt(path))
    data = np.concatenate(data).flatten()
    outdata = np.random.choice(data, number)
    return outdata


alpine_fractions = get_alpine_fractions(return_paths=True)
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
    # to add ashley output I simply need to add the sites here., do this

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


def setup_run_gmp_plus(rch_con, base_mt3d_dir, out_nc, nc_description, ftl_repo=pc5_ftl_repo, dt0=1e4, ttsmax=1e5,
                       ssm_kwargs=None, sft_kwargs=None):
    if ssm_kwargs is None:
        ssm = get_ssm_stress_period_data()
    else:
        ssm = get_ssm_stress_period_data(**ssm_kwargs)

    if sft_kwargs is None:
        sft = get_sft_stress_period_data()
    else:
        sft = get_sft_stress_period_data(**sft_kwargs)

    # get scon from a single run to try and speed up mt3d runs
    scon = None
    # in future I could probably speed up runs by running a single mt3d run.  a test showed a run time of
    # 1412s for a scon of zero and 1335 s for a scon of another run... doesn't seem worthwhile

    setup_run_mt3d_suite(base_mt3d_dir, ftl_repo,
                         ssm_crch={0: rch_con},
                         ssm_stress_period_data={0: ssm},
                         sft_spd={0: sft}, dt0=dt0, ttsmax=ttsmax, scon=scon)

    extract_data(base_mt3d_dir, outfile=out_nc, description=nc_description, nname='mednload', units='g/m3')


def extract_receptor_data(scenario_paths, cbc_paths, outdir):
    """

    :param scenario_paths: dictionary scenario:path to nc file;  scenario must have gmp or cmp in the scenario name
    :param cbc_paths: either the path to the cbc_netcdf (string) or a dictionary scenario and path to teh cbc_netcdf
    :return:
    """
    # extract the raw data from the model runs for all receptors intra and interzone...
    # intrazone data
    intrazone_dir = os.path.join(outdir, 'waimakariri_zone')

    # raw data
    corrected_dir = os.path.join(intrazone_dir, 'corrected_model_data')
    raw_outdir = os.path.join(intrazone_dir, 'raw_model_data')
    if not os.path.exists(corrected_dir):
        os.makedirs(corrected_dir)
    scenarios = scenario_paths.keys()
    pers_names = _np_describe(smt.get_empty_model_grid()).index.values
    outdata = pd.DataFrame(index=sites,
                           columns=pd.MultiIndex.from_product((scenarios, pers_names),
                                                              names=['scenario', 'stat']), dtype=float)
    outdata_raw = pd.DataFrame(index=sites,
                               columns=pd.MultiIndex.from_product((scenarios, pers_names),
                                                                  names=['scenario', 'stat']), dtype=float)
    well_sites = get_well_nums_for_group()
    cmp_waimak_data = {  # the alpine river fractions
        'stream': pd.read_csv(os.path.join(emma_dir, "waimak_per_results_at_points/raw_stocastic_set_str_data.csv"),
                              index_col=0).transpose(),
        'well': pd.read_csv(os.path.join(emma_dir, "waimak_per_results_at_points/raw_stocastic_set_well_data.csv"),
                            index_col=0).transpose()
    }
    for scen, path in scenario_paths.items():
        raw_dir = os.path.join(intrazone_dir, 'raw_model_data', scen)
        if not os.path.exists(raw_dir):
            os.makedirs(raw_dir)
        plot_dir = os.path.join(corrected_dir, '{}_plots'.format(scen))
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        if isinstance(cbc_paths, str):
            use_cbc_path = cbc_paths
        elif isinstance(cbc_paths, dict):
            use_cbc_path = cbc_paths[scen]
        else:
            raise ValueError('unexpected types for cbc_paths: {}'.format(type(cbc_paths)))

        # exctract data from the ucn netcdf and the cbc nettcdf
        get_n_at_points_nc(raw_dir, nsmc_nums=get_stocastic_set(False), ucn_var_name='mednload',
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
                corrected_raw_data = correct_alpine_river(site, cmp_waimak_data, raw_data, well_sites, plot_dir)
            else:
                corrected_raw_data = raw_data.copy()

            # add the stocastic N load
            stocastic_data = add_stocastic_load(site, corrected_raw_data)
            all_raw_data = _np_describe(raw_data)
            for per in pers_names:
                outdata.loc[site, (scen, per)] = stocastic_data[per]
                outdata_raw.loc[site, (scen, per)] = all_raw_data[per]

    outdata.to_csv(os.path.join(corrected_dir, 'all_n_waimak_zone.csv'))
    outdata = outdata.reorder_levels(['stat', 'scenario'], axis=1)
    outdata['1%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_1ths.csv'))
    outdata['5%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_5ths.csv'))
    outdata['10%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_10ths.csv'))
    outdata['25%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_25ths.csv'))
    outdata['50%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_50ths.csv'))
    outdata['75%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_75ths.csv'))
    outdata['95%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_95ths.csv'))
    outdata['99%'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_99ths.csv'))
    outdata['mean'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_mean.csv'))
    outdata['std'].to_csv(os.path.join(corrected_dir, 'n_data_waimak_zone_std.csv'))

    outdata_raw.to_csv(os.path.join(raw_outdir, 'all_n_waimak_zone.csv'))
    outdata_raw = outdata_raw.reorder_levels(['stat', 'scenario'], axis=1)
    outdata_raw['1%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_1ths.csv'))
    outdata_raw['5%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_5ths.csv'))
    outdata_raw['10%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_10ths.csv'))
    outdata_raw['25%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_25ths.csv'))
    outdata_raw['50%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_50ths.csv'))
    outdata_raw['75%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_75ths.csv'))
    outdata_raw['95%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_95ths.csv'))
    outdata_raw['99%'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_99ths.csv'))
    outdata_raw['mean'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_mean.csv'))
    outdata_raw['std'].to_csv(os.path.join(raw_outdir, 'n_data_waimak_zone_std.csv'))

    ##### interzone receptors including n load uncertainty #####
    interzone_outdir = os.path.join(outdir, 'interzone')
    if not os.path.exists(interzone_outdir):
        os.makedirs(interzone_outdir)
    get_interzone_n(scenario_paths,
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


def add_stocastic_load(site, base_data):
    base_data = np.atleast_1d(base_data)
    if site in get_str_ids():
        # stream  # cmp as gmp has a reduction associated
        base_dir = os.path.join(emma_dir, 'nwai_sfeds_nconc_wo_trans_raw_data')
    elif 'wdc' in site:
        # wdc_well # cmp as gmp has a reduction associated
        base_dir = os.path.join(emma_dir,'wdc_well_sfeds_nconc_wo_trans_raw_data')
    else:
        # private well # path must not be ucn as it is more than the maxpath length
        base_dir = os.path.join(emma_dir,'private_well_sfeds_nconc_wo_trans_raw_data')

    modifiers = np.loadtxt(os.path.join(base_dir, '{}.txt'.format(site.replace('wdc_', ''))))

    raw_data = base_data.copy().flatten()
    all_n = raw_data[:, np.newaxis] * modifiers[np.newaxis, :]
    outdata = _np_describe(all_n)
    return outdata


def correct_alpine_river(site, waimak_data, n_data, well_sites, plot_dir):
    print(site)
    warn('correcting data for EMMA results please review plots in: {}'.format(plot_dir))
    if site in get_str_ids():
        # stream stuff
        wai = waimak_data['stream'].loc[site]
        n = n_data
    else:
        wai = waimak_data['well'].loc[well_sites[site]].mean(axis=0)
        n = n_data.mean(axis=0)
        # well stuff
    n.name = 'n'
    wai.name = 'wai'
    plot_data = pd.merge(pd.DataFrame(wai), pd.DataFrame(n),
                         left_index=True, right_index=True)
    model = LR(plot_data.wai, plot_data.n)
    outdata = model.predict(get_alpine_fractions(site))
    model.plot(False, os.path.join(plot_dir, '{}.png'.format(site)))
    return outdata


if __name__ == '__main__':
    for key, paths in alpine_fractions.items():
        for path in paths:
            if not os.path.exists(path):
                print('missing', key, path)

    if False:
        # a test of the extract receptor data... I expect something to break...
        gmp_cbc = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_cbc.nc")
        out_nc_25 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus25per_dairy_ucn.nc")
        out_nc_15 = env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_plus15per_dairy_ucn.nc")
        extract_receptor_data(scenario_paths={
            'gmp_15_red': out_nc_15,
            'gmp_25_red': out_nc_25,
        },
            cbc_paths={
                'gmp_15_red': gmp_cbc,
                'gmp_25_red': gmp_cbc,
            },
            outdir=r"C:\Users\matth\Downloads\test_extract_receptor_data")

    print('done')
