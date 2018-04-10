# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 10/04/2018 8:45 AM
"""

from __future__ import division
from core import env
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import os
import pandas as pd


def get_current_pathway_n(mode):
    if mode == '50th':
        private_key = '50%_CP_con'
        wdc_key = '50%_gmp_con'
        stream_key = 'Median Current Pathways'
    elif mode == '95th':
        private_key = '95%_CP_con'
        wdc_key = '95%_gmp_con'
        stream_key = '95th percentile Current Pathways'
    else:
        raise ValueError('unexpected value for mode: {} expected "50th" or "95th"'.format(mode))
    outdata = {}

    # wells
    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N '
                         r'results wells summary.xlsx',
                         sheetname='CP WDC', index_col=0)
    for key in wdc_wells:
        outdata[key] = data.loc[key, wdc_key]

    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N '
                         r'results wells summary.xlsx',
                         sheetname='CP private', index_col=0)
    for key in private_wells:
        outdata[key] = data.loc[key, private_key]

    # steams
    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model '
                         r'simulations and results\ex_bd_va\n_results\N results for northern tribs ZC workshop\N load '
                         r'reductions Streams 050418.xlsx', skiprows=1, index_col=0)

    for key in streams:
        outdata[key] = data.loc[str_header_conversion[key], stream_key].iloc[0]

    return outdata


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


def gen_well_targets(scenario):
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
    elif scenario == 'permissive':
        for key in wdc_wells:
            outdata[key] = 5.65
        for key in private_wells:
            outdata[key] = 7.1  # a number which should represent close to below the drinking water standard
    elif scenario == 'current_measured':
        for key in private_wells:
            outdata[key] = 3.5
            outdata.update({
                'wdc_Cust': 4.8,  # estimate from nearby, no measured wdc  n
                'wdc_Fernside': 1.1,
                'wdc_Kaiapoi': 1.5,
                'wdc_Kairaki': 0.4,
                'wdc_Mandeville': 2.1,
                'wdc_Ohoka': 3.0,  # estimate from nearby, no measured wdc  n
                'wdc_Oxford Urban': 2.7,
                'wdc_Pegasus': 0.05,
                'wdc_Poyntzs Road': 3.0,  # estimate from nearby, no measured wdc  n
                'wdc_Rangiora': 1.5,
                'wdc_Waikuku': 0.7,
                'wdc_West Eyreton': 3.0,  # estimate from nearby, no measured wdc n
            })
    elif scenario == 'zonal_average':
        for key in (wdc_wells | private_wells):
            outdata[key] = 3.56
    else:
        raise NotImplementedError('scenario: {} has not been implemented'.format(scenario))

    return outdata


def gen_stream_targets(scenario):
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

    return outdata


def calc_per_reduction_rasters(outdir, name, mode, well_targets, stream_targets, waimak_target=27,
                               interzone_target=None, include_interzone=False):
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
    :return: save rasters percentage reduction needed (0-100) and limiting raster:
                             {1:'interzone', 2:'wdc_wells', 3:'private_wells', 4: 'springfed_streams', 5: 'waimakariri'
                             6: 'all same'}
    """
    # something to pull togeather different well/stream target scenarios and look at the percentage reduction required for
    # each
    # output as rasters
    # also a map of limiting factors
    # include waimak target of 0.1 as given.
    # todo is waimak a 27% reduction from the 95th assuming here this is constant
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    current_paths = get_current_pathway_n(mode)
    base_shp_path = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results"
                            r"\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons")

    # interzone
    if include_interzone:
        raise NotImplementedError
    else:
        interzone_reduction = smt.get_empty_model_grid() * np.nan

    # wdc wells
    wdc_dir = os.path.join(base_shp_path, 'wdc_wells')
    wdc_reductions = []
    for well in wdc_wells:
        temp = smt.shape_file_to_model_array(os.path.join(wdc_dir, well.replace('wdc_', '') + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - well_targets[well] / current_paths[well]) * 100
        wdc_reductions.append(temp[np.newaxis])
    wdc_reductions = np.nanmax(np.concatenate(wdc_reductions, axis=0), axis=0)

    # private wells
    private_dir = os.path.join(base_shp_path, 'private_wells')
    private_reductions = []
    for well in private_wells:
        temp = smt.shape_file_to_model_array(os.path.join(private_dir, well + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - well_targets[well] / current_paths[well]) * 100
        private_reductions.append(temp[np.newaxis])
    private_reductions = np.nanmax(np.concatenate(private_reductions, axis=0), axis=0)

    # streams
    stream_dir = os.path.join(base_shp_path, 'second_tranche')
    stream_reductions = []
    for stream in streams:
        temp = smt.shape_file_to_model_array(os.path.join(stream_dir, stream + '.shp'),
                                             'Id', True)
        temp[np.isfinite(temp)] = (1 - stream_targets[stream] / current_paths[stream]) * 100
        stream_reductions.append(temp[np.newaxis])
    stream_reductions = np.nanmax(np.concatenate(stream_reductions, axis=0), axis=0)

    # waimakariri
    waimak_reduction = smt.shape_file_to_model_array(os.path.join(base_shp_path, r"waimakariri_river\waimakariri.shp"),
                                                     'Id', True)
    waimak_reduction[np.isfinite(waimak_reduction)] = waimak_target

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
    reason = np.nanargmax(temp.astype(int), axis=0).astype(
        float) + 1
    reason[idx & np.isfinite(reduction)] = 6
    reason[np.isnan(reduction)] = np.nan
    smt.array_to_raster(os.path.join(outdir, '{}_reason.tif'.format(name)), reason)
    smt.array_to_raster(os.path.join(outdir, '{}_reduction.tif'.format(name)), reduction)


if __name__ == '__main__':
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
                                   interzone_target=None, include_interzone=False)
