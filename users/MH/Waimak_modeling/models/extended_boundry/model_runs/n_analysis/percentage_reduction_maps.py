# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 10/04/2018 8:45 AM

original percentage reduction maps.  scripts have been modified to include pc5pa rules, some more scenarios, and mar
options, but the main component has not been modified.
"""

from __future__ import division
from core import env
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import os
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.n_load_layers import \
    get_gmp_load_raster, get_pc5pa_additonal_load


def get_current_pathway_n(mode, conservative_zones):
    if not isinstance(mode, dict):
        if mode == '50th':
            mode = {}
            for key in private_wells:
                mode[key] = '50%_gmp_con'
            for key in wdc_wells:
                mode[key] = '50%_gmp_con'
            for key in streams:
                mode[key] = 'Median Current Pathways'
        elif mode == '95th':
            mode = {}
            for key in private_wells:
                mode[key] = '95%_gmp_con'
            for key in wdc_wells:
                mode[key] = '95%_gmp_con'
            for key in streams:
                mode[key] = '95th percentile Current Pathways'
        else:
            raise ValueError('unexpected value for mode: {} expected "50th" or "95th"'.format(mode))
    outdata = {}

    # wells
    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N '
                         r'results wells summary.xlsx',
                         sheetname='GMP WDC', index_col=0)
    for key in wdc_wells:
        outdata[key] = data.loc[key, mode[key]]

    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N '
                         r'results wells summary.xlsx',
                         sheetname='GMP private', index_col=0)
    for key in private_wells:
        outdata[key] = data.loc[key, mode[key]]

    # steams
    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N load '
                         r'reductions Streams 050418.xlsx', skiprows=1, index_col=0)

    for key in streams:
        outdata[key] = data.loc[str_header_conversion[key], mode[key]].iloc[0]

    # add PA n load increases
    pa_n = get_pa_reductions(conservative_zones)
    for key in set(outdata.keys()) - streams:  # streams are not included as they don't change size
        outdata[key] += outdata[key] * pa_n[key] / 100

    return outdata


def get_pa_reductions(conservative_zones):
    # get teh pa reductions to apply to the below script. need to integrate the conservative and non conservitive zones
    # pull from my PA load scripts
    outdata = {}

    # streams
    stream_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                          r"s\ex_bd_va\n_results\nwaimak_springfeds\pa_rules_nwaimak_springfeds\load_overviews_"
                          r"with_paN.csv")
    data = pd.read_csv(stream_path, index_col=0)
    data.loc[:, 'pa_increase'] = data.loc[:, 'total_pa_N_kg'] / data.loc[:, 'gmp_nload_kg'] * 100
    for key in streams:
        outdata[key] = data.loc[key, 'pa_increase']

    # wdc wells
    if conservative_zones == 'conservative':
        wdc_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd"
                           r"_va\n_results\wdc_wells\pa_rules_wdc_wells\load_overviews_with_paN.csv")
    elif conservative_zones == 'permissive':
        wdc_path = env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_"
            r"va\n_results\wdc_wells_90\pa_rules_wdc_wells_90\load_overviews_with_paN.csv")
    elif conservative_zones == 'use_mix':
        wdc_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex"
                           r"_bd_va\n_results\wdc_use_mix\pa_rules_wdc_use_mix\load_overviews_with_paN.csv")
    else:
        raise NotImplementedError('conservitive zones not implemented: {}'.format(conservative_zones))

    data = pd.read_csv(wdc_path, index_col=0)
    data.loc[:, 'pa_increase'] = data.loc[:, 'total_pa_N_kg'] / data.loc[:, 'gmp_nload_kg'] * 100
    for key in wdc_wells:
        outdata[key] = data.loc[key.replace('wdc_', ''), 'pa_increase']

    # private wells
    if conservative_zones == 'conservative':
        private_path = env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_"
            r"bd_va\n_results\private_wells\pa_rules_private_wells\load_overviews_with_paN.csv")
    elif conservative_zones == 'permissive':
        private_path = env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex"
            r"_bd_va\n_results\private_wells_90\pa_rules_private_wells_90\load_overviews_with_paN.csv")
    elif conservative_zones == 'use_mix':
        private_path = env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex"
            r"_bd_va\n_results\private_wells_90\pa_rules_private_wells_90\load_overviews_with_paN.csv")
    else:
        raise NotImplementedError('conservitive zones not implemented: {}'.format(conservative_zones))

    data = pd.read_csv(private_path, index_col=0)
    data.loc[:, 'pa_increase'] = data.loc[:, 'total_pa_N_kg'] / data.loc[:, 'gmp_nload_kg'] * 100
    for key in private_wells:
        outdata[key] = data.loc[key, 'pa_increase']

    for key in outdata.keys():
        outdata[key] *= 0.5  # current pathways estmates a rough uptake of 50%
    return outdata


stream_apply_mar = {'cam_bramleys_s': False,
                    'courtenay_kaiapoi_s': True,
                    'cust_skewbridge': True,
                    'kaiapoi_harpers_s': True,
                    'kaiapoi_island_s': True,
                    'ohoka_island_s': True
                    }

wdc_current_measured = {
    'wdc_Cust': 0.3,
    'wdc_Fernside': 1.1,
    'wdc_Kaiapoi': 1.5,
    'wdc_Kairaki': 0.4,
    'wdc_Mandeville': 3.0,
    'wdc_Ohoka': 8.0,
    'wdc_Oxford Urban': 2.7,
    'wdc_Pegasus': 0.02,
    'wdc_Poyntzs Road': 9.0,
    'wdc_Rangiora': 1.5,
    'wdc_Waikuku': 0.7,
    'wdc_West Eyreton': 2.0,
}

wdc_wells = {'wdc_Cust',
             'wdc_Fernside',
             'wdc_Kaiapoi',
             'wdc_Kairaki',
             'wdc_Mandeville',
             'wdc_Ohoka',
             'wdc_Oxford Urban',
             'wdc_Pegasus',
             'wdc_Poyntzs Road',
             'wdc_Rangiora',
             'wdc_Waikuku',
             'wdc_West Eyreton'}

private_wells = {'Clarkville',
                 'Cust',
                 'Eyreton_deep',
                 'Eyreton_shallow',
                 'Fernside',
                 'Flaxton',
                 'Horellville',
                 'Mandeville',
                 'North East Eyrewell_deep',
                 'North East Eyrewell_shallow',
                 'North West Eyrewell_deep',
                 'North West Eyrewell_shallow',
                 'Ohoka_deep',
                 'Ohoka_shallow',
                 'Rangiora',
                 'Springbank',
                 'Summerhill',
                 'Swannanoa_deep',
                 'Swannanoa_shallow',
                 'Waikuku',
                 'West Eyreton_deep',
                 'West Eyreton_shallow',
                 'Woodend - Tuahiwi'}

streams = {'cam_bramleys_s',
           'courtenay_kaiapoi_s',
           'cust_skewbridge',
           'kaiapoi_harpers_s',
           'kaiapoi_island_s',
           'ohoka_island_s'}

str_header_conversion = {'cam_bramleys_s': 'Cam  at Bramleys Rd',
                         'courtenay_kaiapoi_s': 'Courtenay  at Kaiapoi confluence',
                         'cust_skewbridge': 'Cust  at Skewbridge Rd',
                         'kaiapoi_harpers_s': 'Kaiapoi  at Harpers Rd',
                         'kaiapoi_island_s': 'Kaiapoi  at Island Rd',
                         'ohoka_island_s': 'Ohoka  at Island  at Rd'}


def get_mode(scenario):
    if scenario == 'least_pain':
        mode = '50th'
    elif scenario == 'middle_option':
        mode = {}
        for key in (wdc_wells | private_wells):
            mode[key] = '50%_gmp_con'
        for key in streams:
            mode[key] = '95th percentile Current Pathways'
        mode['cam_bramleys_s'] = 'Median Current Pathways'
    elif scenario == 'most_gain':
        mode = '95th'
    else:
        raise NotImplementedError('scenario: {} not implmented'.format(scenario))

    return mode


def gen_waimak_targets(scenario):
    if scenario == 'waimak_None':
        out = 0
    elif scenario == 'least_pain':
        out = 0
    elif scenario == 'middle_option':
        out = 28
    elif scenario == 'most_gain':
        out = 28
    else:
        raise NotImplementedError('scenario not implemented: {}'.format(scenario))
    return out


def gen_well_targets(scenario, wdc_none=False, private_none=False):
    """
    get the well concentration targets for different scenarios
    :param scenario: a keyword for different scenarios
                        'half_mav': set all well concentrations to 5.65 mg/l
                        'current_measured': set well concentrations to current measured
                                            or zonal average of 3.5 mg/l for private supply wells
                        'zonal_average': set all well concentrations to 3.5 mg/l
    :return:
    """
    # some helper script to rapidly handle the private and wdc wells.
    # half MAV
    # current 'measured'
    outdata = {}
    if scenario == 'half_mav':
        for key in (wdc_wells | private_wells):
            outdata[key] = 5.65

    # following 3 scenarios match zeb's options
    elif scenario == 'least_pain':
        for key in wdc_wells:
            outdata[key] = 5.65
        for key in private_wells:
            outdata[key] = 7.1  # a number which should represent close to below the drinking water standard
    elif scenario == 'middle_option':
        for key in wdc_wells:
            outdata[key] = 5.65
        for key in private_wells:
            outdata[key] = 5.65  # a number which should represent close to below the drinking water standard
    elif scenario == 'most_gain':
        for key in wdc_wells:
            outdata[key] = 5.65
        for key in private_wells:
            outdata[key] = 5.65  # a number which should represent close to below the drinking water standard

    elif scenario == 'permissive':
        for key in wdc_wells:
            outdata[key] = 5.65
        for key in private_wells:
            outdata[key] = 7.1  # a number which should represent close to below the drinking water standard
    elif scenario == 'current_measured':
        for key in private_wells:
            outdata[key] = 3.5
            outdata.update(wdc_current_measured)
    elif scenario == 'zonal_average':
        for key in (wdc_wells | private_wells):
            outdata[key] = 3.56
    else:
        raise NotImplementedError('scenario: {} has not been implemented'.format(scenario))

    if wdc_none:  # to allow not passing targets for private wells
        for key in (wdc_wells):
            outdata[key] = 1e30
    if private_none:
        for key in (private_wells):
            outdata[key] = 1e30

    return outdata


def gen_stream_targets(scenario, stream_none=False):
    """
    # some helper function to sort out the stream targets
    :param scenario: a keyword for different scenarios
                     'cultural': 1 mg/l target for all bodies
                     'current_measured': the targets at the current measured values
                     'top_current_nof': the maximum allowable under the NPS (e.g. the top of the current
                     measured NOF band or the national bottom line of 6.9)
                     'reasonable_reduction': the reasonable reduction to the values of the NOF bands,
                                             1, 2.4, 3.8(half cband), 6.9 more details in the implementation
    :return:
    """
    outdata = {}
    if scenario == 'cultural':
        for key in streams:
            outdata[key] = 1
    # following 3 scenarios match zeb's options
    elif scenario == 'least_pain':
        outdata.update({'cam_bramleys_s': 2.4,
                        'courtenay_kaiapoi_s': 6.9,
                        'cust_skewbridge': 6.9,
                        'kaiapoi_harpers_s': 6.9,  # current measured is 9.4, set to national bottom line
                        'kaiapoi_island_s': 6.9,
                        'ohoka_island_s': 6.9})
    elif scenario == 'middle_option':
        outdata.update({'cam_bramleys_s': 1.6,
                        'courtenay_kaiapoi_s': 6.9,
                        'cust_skewbridge': 4.7,
                        'kaiapoi_harpers_s': 6.9,  # current measured is 9.4, set to national bottom line
                        'kaiapoi_island_s': 6.9,
                        'ohoka_island_s': 6.9})
    elif scenario == 'most_gain':
        outdata.update({'cam_bramleys_s': 1.5,
                        'courtenay_kaiapoi_s': 3.1,
                        'cust_skewbridge': 3.8,
                        'kaiapoi_harpers_s': 3.8,  # current measured is 9.4, set to national bottom line
                        'kaiapoi_island_s': 5.4,
                        'ohoka_island_s': 4.5})


    elif scenario == 'current_measured':
        outdata.update({'cam_bramleys_s': 1.5,
                        'courtenay_kaiapoi_s': 2.5,
                        'cust_skewbridge': 4.7,
                        'kaiapoi_harpers_s': 6.9,  # current measured is 9.4, set to national bottom line
                        'kaiapoi_island_s': 5.4,
                        'ohoka_island_s': 4.5})

    elif scenario == 'top_current_nof':
        outdata.update({'cam_bramleys_s': 2.4,
                        'courtenay_kaiapoi_s': 6.9,  # poor measurments
                        'cust_skewbridge': 6.9,
                        'kaiapoi_harpers_s': 6.9,
                        'kaiapoi_island_s': 6.9,
                        'ohoka_island_s': 6.9})

    elif scenario == 'asperational_reduction':
        outdata.update({'cam_bramleys_s': 1,  # reduction to a band
                        'courtenay_kaiapoi_s': 3.8,  # aim to top of B band
                        'cust_skewbridge': 3.8,  # aim for half C band
                        'kaiapoi_harpers_s': 3.8,  # aim for national bottom line
                        'kaiapoi_island_s': 3.8,  # aim for half C band
                        'ohoka_island_s': 3.8})  # aim for half C band

    elif scenario == 'attainable_reduction':
        outdata.update({'cam_bramleys_s': 1,  # reduction to a band
                        'courtenay_kaiapoi_s': 6.9,  # aim to top of c band
                        'cust_skewbridge': 3.8,  # aim for half C band
                        'kaiapoi_harpers_s': 6.9,  # aim for national bottom line
                        'kaiapoi_island_s': 6.9,  # aim for top C band
                        'ohoka_island_s': 6.9})  # aim for top C band

    else:
        raise NotImplementedError('scenario: {} is not implemented'.format(scenario))

    if stream_none:  # to allow passing no target for streams
        for key in streams:
            outdata[key] = 1e30

    return outdata


def get_interzone_reduction(target_load=8, pc5pa=False):
    """
    modelled off GMP loads assumes the conservative zone
    :param target_load: the load that we want to achieve
    :return:
    """

    gmp_load = get_gmp_load_raster()

    if pc5pa:
        gmp_load += get_pc5pa_additonal_load()

    interzone_reductions = target_load / gmp_load
    shp_file_path = env.sci("Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and resul"
                            "ts\ex_bd_va\capture_zones_particle_tracking\source_zone_polygon"
                            "s\interzone\conservative_interzone.shp")
    interzone_idx = np.isfinite(smt.shape_file_to_model_array(shp_file_path, 'Id', True))
    interzone_reductions[~interzone_idx] = 0
    return interzone_reductions


def calc_per_reduction_rasters(outdir, name, mode, well_targets, stream_targets, waimak_target=27,
                               mar_percentage=0, pc5_pa_rules=False, conservative_shp='conservative',
                               interzone_target_load=None, include_interzone=False, save_reason=True):
    """
    pull togeather the different current pathway results at 95th percentile and calculate reduction rasters (from 95th),
    also calculate a 5 integer raster that says whether it is surface water or ground water (split into wdc and private)
    concentration or equal limited or interzone

    :param outdir: directory to save rasters
    :param name: name to identify the pathways
    :param mode: which current pathways to use '50th','95th'
    :param well_targets: the well targets as a dictionary
    :param stream_targets: dictionary for stream targets
    :param waimak_target: percentage reduction beyond gmp default = 27%
    :param mar_percentage: percentage handled by mar
    :param pc5_pa_rules:  apply reductions based on setting pc5pa rules to (0,0)
    :param conservative_shp: boolean if True, use the larger zones if False use the smaller zones
    :param save_reason: boolean if True will save the reason, otherwise will not
    :return: save rasters percentage reduction needed (0-100) and limiting raster:
                             {1:'interzone', 2:'wdc_wells', 3:'private_wells', 4: 'springfed_streams', 5: 'waimakariri'
                             6: 'all same'}
    """
    # something to pull togeather different well/stream target scenarios and look at the percentage reduction required for
    # each
    # output as rasters
    # also a map of limiting factors
    # include waimak target of 0.1 as given.
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    current_paths = get_current_pathway_n(mode, conservative_shp)
    base_shp_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results"
                            r"\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons")

    if pc5_pa_rules:
        pa_reductions = get_pa_reductions(conservative_shp)
    # interzone
    if include_interzone:
        interzone_reduction = get_interzone_reduction(interzone_target_load, pc5_pa_rules)
    else:
        interzone_reduction = smt.get_empty_model_grid() * np.nan

    # wdc wells
    if conservative_shp == 'conservative':
        wdc_dir = os.path.join(base_shp_path, 'wdc_wells')
    elif conservative_shp == 'permissive':
        wdc_dir = os.path.join(base_shp_path, 'wdc_wells_90_named_right')
    elif conservative_shp == 'use_mix':
        wdc_dir = os.path.join(base_shp_path, 'wdc_use_mix')  # a mix of the 90 for most zones unless serve > 1000 ppl
    else:
        raise NotImplementedError('conservitive shape not implemented: {}'.format(conservative_shp))

    wdc_reductions = []
    for well in wdc_wells:
        temp = smt.shape_file_to_model_array(os.path.join(wdc_dir, well.replace('wdc_', '') + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - well_targets[well] / current_paths[well]) * 100
        if pc5_pa_rules:
            temp[np.isfinite(temp)] += - pa_reductions[well]

        wdc_reductions.append(temp[np.newaxis])
    wdc_reductions = np.nanmax(np.concatenate(wdc_reductions, axis=0), axis=0)

    # private wells
    if conservative_shp == 'conservative':
        private_dir = os.path.join(base_shp_path, 'private_wells')
    elif conservative_shp == 'permissive':
        private_dir = os.path.join(base_shp_path, 'private_wells_90_named_right')
    elif conservative_shp == 'use_mix':
        private_dir = os.path.join(base_shp_path, 'private_wells_90_named_right')
    else:
        raise NotImplementedError('conservitive shape not implemented: {}'.format(conservative_shp))

    private_reductions = []
    for well in private_wells:
        temp = smt.shape_file_to_model_array(os.path.join(private_dir, well + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - well_targets[well] / current_paths[well]) * 100
        if pc5_pa_rules:
            temp[np.isfinite(temp)] += - pa_reductions[well]
        private_reductions.append(temp[np.newaxis])
    private_reductions = np.nanmax(np.concatenate(private_reductions, axis=0), axis=0)

    # streams
    stream_dir = os.path.join(base_shp_path, 'second_tranche')
    stream_reductions = []
    for stream in streams:
        temp = smt.shape_file_to_model_array(os.path.join(stream_dir, stream + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - stream_targets[stream] / current_paths[stream]) * 100 - (mar_percentage *
                                                                                                stream_apply_mar[
                                                                                                    stream])
        if pc5_pa_rules:
            temp[np.isfinite(temp)] += - pa_reductions[well]
        stream_reductions.append(temp[np.newaxis])
    stream_reductions = np.nanmax(np.concatenate(stream_reductions, axis=0), axis=0)

    # waimakariri
    waimak_reduction = smt.shape_file_to_model_array(os.path.join(base_shp_path, r"waimakariri_river\waimakariri.shp"),
                                                     'Id', True)
    waimak_reduction[np.isfinite(waimak_reduction)] = waimak_target
    if pc5_pa_rules:
        waimak_reduction += -1

    temp = np.concatenate([interzone_reduction[np.newaxis],
                           wdc_reductions[np.newaxis],
                           private_reductions[np.newaxis],
                           stream_reductions[np.newaxis],
                           waimak_reduction[np.newaxis]])
    idx = np.all([e.astype(int) == interzone_reduction.astype(int) for e in [wdc_reductions,
                                                                             private_reductions,
                                                                             stream_reductions,
                                                                             waimak_reduction]], axis=0)
    reduction = np.nanmax(temp, axis=0)
    reduction[reduction < 0] = 0
    reason = np.nanargmax(temp.astype(int), axis=0).astype(
        float) + 1
    reason[idx & np.isfinite(reduction)] = 6
    reason[np.isnan(reduction)] = np.nan
    if save_reason:
        smt.array_to_raster(os.path.join(outdir, '{}_reason.tif'.format(name)), reason)
    smt.array_to_raster(os.path.join(outdir, '{}_reduction.tif'.format(name)), reduction)


if __name__ == '__main__':
    # not changed since original (other than this line)
    import itertools

    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_reductions_excl_interzone"

    for well, stream, mode in itertools.product(['half_mav', 'permissive', 'current_measured'],
                                                ['cultural', 'current_measured', 'top_current_nof',
                                                 'asperational_reduction', 'attainable_reduction'],
                                                ['50th', '95th']):
        print((well, stream, mode))
        well_targets = gen_well_targets(well)
        stream_targets = gen_stream_targets(stream)

        calc_per_reduction_rasters(outdir,
                                   'cp_{}_well_{}_stream_{}'.format(mode, well, stream),
                                   mode, well_targets, stream_targets, waimak_target=27,
                                   interzone_target_load=None, include_interzone=False)
