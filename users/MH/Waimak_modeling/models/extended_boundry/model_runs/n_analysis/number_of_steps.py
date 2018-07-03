# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/06/2018 11:39 AM
"""

from __future__ import division
from core import env
import geopandas as gpd
import pandas as pd
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.run_source_uncertainty import \
    spatial_overlays
from percentage_reduction_maps import gen_stream_targets, gen_well_targets, get_current_pathway_n, gen_interzone_targets
import pickle
import os
from copy import deepcopy
import itertools
import numpy as np

use_soil_mapper = {
    'dairy_lowdrain': 1,
    'sbda_lowdrain': 1,
    'dairy_highdrain': 0,
    'sbda_highdrain': 0,
}
use_luse_mapper = {
    'dairy_lowdrain': 1,
    'sbda_lowdrain': 2,
    'dairy_highdrain': 1,
    'sbda_highdrain': 2,
}

order = ['Streams', 'cam_bramleys_s', 'courtenay_kaiapoi_s', 'kaiapoi_harpers_s', 'kaiapoi_island_s',
         'ohoka_island_s', 'Private Wells', 'Clarkville', 'Cust', 'cust_skewbridge', 'Eyreton_deep',
         'Eyreton_shallow', 'Fernside', 'Flaxton', 'Horellville', 'Mandeville', 'North East Eyrewell_deep',
         'North East Eyrewell_shallow', 'North West Eyrewell_deep', 'North West Eyrewell_shallow', 'Ohoka_deep',
         'Ohoka_shallow', 'Rangiora', 'Springbank', 'Summerhill', 'Swannanoa_deep', 'Swannanoa_shallow', 'Waikuku',
         'West Eyreton_deep', 'West Eyreton_shallow', 'Woodend - Tuahiwi', 'Wdc Wells', 'wdc_Cust', 'wdc_Fernside',
         'wdc_Kaiapoi', 'wdc_Kairaki', 'wdc_Mandeville', 'wdc_Ohoka', 'wdc_Oxford Urban', 'wdc_Pegasus',
         'wdc_Poyntzs Road', 'wdc_Rangiora', 'wdc_Waikuku', 'wdc_West Eyreton', 'Interzone', 'conservative_interzone',
         'highly_likely_interzone']


lclasses = {'dairy_lowdrain',
            'sbda_lowdrain',
            'dairy_highdrain',
            'sbda_highdrain',
            'lifestyle'}


def gen_site_shape_dict():
    sites = gen_well_targets('preferred')
    sites.update(gen_stream_targets('preferred'))
    sites = sites.keys()
    outdata = {}
    paths = glob(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd"
                 r"_va\capture_zones_particle_tracking\source_zone_"
                 r"polygons\private_wells_90_named_right\*.shp")
    paths.extend(glob(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd"
                      r"_va\capture_zones_particle_tracking\source_zone_"
                      r"polygons\wdc_use_mix\*.shp"))
    paths.extend(glob(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd"
                      r"_va\capture_zones_particle_tracking\source_zone_"
                      r"polygons\second_tranche\*.shp"))
    paths = pd.Series(paths)

    for site in sites:
        if 'wdc' in site:
            idx = paths.str.contains(site.replace('wdc_', '')) & paths.str.contains('wdc_use_mix')
        else:
            idx = paths.str.contains(site) & ~paths.str.contains('wdc_use_mix')

        outdata[site] = paths.loc[idx].iloc[0]

    outdata['conservative_interzone'] = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\interzone\conservative_interzone.shp"
    outdata['highly_likely_interzone'] = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\interzone\highly_likely_interzone.shp"

    return outdata


def get_average_n(recalc=False):
    pickle_path = r"T:\Temp\temp_gw_files\average_n.p"
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    base_nload_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                              'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')

    land_types = {'Arable': 2,  # mean of 18 kg/ha/year
                  'DairyFarm': 1,  # mean 36 kg/ha/year
                  'DairySupport': 1,  # mean of 24 kg/ha/year
                  'ForestTussock': 0,  # mean of 1 kg/ha/year
                  'Horticulture': 2,  # mean 14 kg/ha/year
                  'Lifestyle': 3,  # mean of 5kg/ha/year
                  'NotFarm': 0,  # all 2 kg/ha/year
                  'OtherinclGolf': 0,  # mean of 11.8 kg/ha/year
                  'Pigs': 2,  # mean of 28kg/ha/year
                  'SheepBeefDeer': 2,  # mean of 12.5 kg/ha/year
                  'SheepBeefHill': 2,  # mean of 3kg/ha/year
                  'Unknown': 2}  # mean of 18.49 kg/ha/year

    soil_classes = {'': 0,
                    'WW': 0,
                    'L': 0,
                    'XL': 0,
                    'D': 0,
                    'Pd': 1,
                    'PdL': 1,
                    'F1': 0,
                    'S3': 0,
                    'M': 0,
                    'VL': 0,
                    'S1': 0,
                    'F3': 0,
                    'S2': 0,
                    'S4': 0,
                    'F2': 0,
                    'O': 0
                    }

    base_data = gpd.read_file(base_nload_path)
    base_data.replace({'luscen_cat': {
        'Forest-Tussock': 'ForestTussock',
        'Other-inclGolf': 'OtherinclGolf',
        'Sheep-Beef-Deer': 'SheepBeefDeer',
        'SheepBeef-Hill': 'SheepBeefHill',
    }}, inplace=True)
    base_data.loc[base_data.MGMSoilCod.isnull(), 'MGMSoilCod'] = ''

    base_data.loc[:, 'use_soil'] = base_data.MGMSoilCod.replace(soil_classes)
    base_data.loc[:, 'use_ltype'] = base_data.luscen_cat.replace(land_types)

    site_paths = gen_site_shape_dict()
    outdata = {}
    for site, path in site_paths.items():
        temp = gpd.read_file(path)
        temp_n = spatial_overlays(base_data, temp)
        temp_n.loc[:, 'use_area'] = temp_n.geometry.area
        temp_out = {}
        for lclass in lclasses:
            if lclass == 'lifestyle':
                idx = temp_n.use_ltype == 3
            else:
                idx = (temp_n.use_ltype == use_luse_mapper[lclass]) & (temp_n.use_soil == use_soil_mapper[lclass])
            total_area = temp_n.loc[idx, 'use_area'].sum()
            if total_area == 0:
                n = 0
            else:
                n = (temp_n.loc[idx, 'use_area'] * temp_n.loc[idx, 'nconc_gmp']).sum() / total_area
            temp_out[lclass] = n
        outdata[site] = temp_out
    pickle.dump(outdata, open(pickle_path, 'w'))
    return outdata


def get_fractions():
    # this comes from a series of mt3d runs
    base_path = env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\landuse_drainage_fractions\waimakariri_zone\corrected_model_data\n_data_waimak_zone_50ths.csv")

    basedata = pd.read_csv(base_path, index_col=0)
    interzone_data = pd.read_csv(r'P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and resu'
                                 r'lts\ex_bd_va\zc_n_sols\landuse_drainage_fractions\interzone\n_'
                                 r'data_interzone_50ths.csv',
                                 index_col=[0,1])
    basedata.loc['conservative_interzone'] = interzone_data.loc[('deep unnamed', 'full_city')]
    basedata.loc['highly_likely_interzone'] = interzone_data.loc[('deep unnamed', 'full_city')]
    basedata.loc[:, 'other'] = 1 - basedata.sum(axis=1)
    return basedata.transpose().to_dict()


def number_of_steps(step_reductions, pa_00, target_scheme):
    """
    calculate the number of steps to reach a selected target scheme
    :param step_reductions: a dictionary of the percentage reduction to apply at each step e.g.
                            {'dairy_lowdrain' :20,
                            'sbda_lowdrain':10,
                            'dairy_highdrain':20,
                            'sbda_highdrain':5,
                            'lifestyle':0}
    :param pa_00: boolean if True, set teh pa rules to 0,0
    :param target_scheme: one of: ['preferred', 'alternate']
    :return:
    """
    assert isinstance(step_reductions, dict)
    assert set(step_reductions.keys()) == lclasses
    use_step_reductions = {}

    for k, v in step_reductions.items():
        use_step_reductions[k] = 1 - v / 100

    targets = gen_well_targets(target_scheme)
    targets.update(gen_stream_targets(target_scheme))
    targets.update(gen_interzone_targets(target_scheme))
    mode = {k: '50%' for k in targets.keys()}
    current_paths = get_current_pathway_n(mode=mode, conservative_zones='use_mix', from_mt3d_runs=True,
                                          mt3d_add_pa=not pa_00, inc_interzone=True)

    average_n = get_average_n()
    fractions = get_fractions()
    outdata = pd.DataFrame(index=targets.keys())
    for site in targets.keys():
        # calculate the average n concentration for the missing components (check this)
        n_proj = deepcopy(current_paths[site])
        if n_proj < targets[site]:
            continue
        n_proj_temp = deepcopy(current_paths[site])
        for lclass in lclasses:
            n_proj_temp += - average_n[site][lclass] * fractions[site][lclass]

        other_n_con = n_proj_temp / fractions[site]['other']
        outdata.loc[site, 'other_n_conc'] = other_n_con

        # apply reduction, calculate new concentration
        new_n = other_n_con * fractions[site]['other']
        for lclass in lclasses:
            new_n += average_n[site][lclass] * use_step_reductions[lclass] * fractions[site][lclass]
        outdata.loc[site, 'step_1_con'] = new_n

        # from new concentration calculate the number of steps to the target
        dif = n_proj - new_n
        num_steps = (n_proj - targets[site]) / dif
        if num_steps * max(step_reductions.values()) > 100:
            num_steps = 999 # more steps than possible

        outdata.loc[site, 'nsteps'] = num_steps

    return outdata

def run_scenarios(outdir): #todo run scenarios and make timelines # set max threshold to 100%
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdata = {}
    target_options = ['preferred', 'alternate']
    pas = [True, False]

    # make reduction options
    dairy_options = [10, 25]
    sheep_options = [0, 5]
    lifestyle_options = [0,2]
    include_ld = [True, False]

    reduction_options = {}
    for d, s, ly, ld in itertools.product(dairy_options, sheep_options, lifestyle_options, include_ld):
        ld_name = 'exc'
        if ld:
            ld_name = 'inc'

        name = 'd{}_s{}_L{}_{}ld'.format(d, s, ly, ld_name)
        temp = {}
        temp['dairy_highdrain'] = d
        temp['sbda_highdrain'] = s
        temp['lifestyle'] = ly

        if ld:
            temp['dairy_lowdrain'] = d
            temp['sbda_lowdrain'] = s
        else:
            temp['dairy_lowdrain'] = 0
            temp['sbda_lowdrain'] = 0

        reduction_options[name] = temp

    for red, pc5pa00, tar in itertools.product(reduction_options.keys(), pas, target_options):
        pa_name = 'without'
        if pc5pa00:
            pa_name = 'with'
        outname = '{}_{}_pc5pa00_{}_tar'.format(red, pa_name, tar)

        temp_data = number_of_steps(step_reductions=reduction_options[red],
                                    pa_00=pc5pa00,
                                    target_scheme=tar)
        temp_data.to_csv(os.path.join(outdir, '{}.csv'.format(outname)))
        outdata[outname] = temp_data.nsteps.round(1)
    outdata = pd.DataFrame(outdata)
    outdata.loc['Streams', :] = np.nan
    outdata.loc['Wdc Wells', :] = np.nan
    outdata.loc['Private Wells', :] = np.nan
    outdata = outdata.loc[order]
    outdata.to_csv(os.path.join(outdir, 'all_nsteps_summary.csv'))
    return outdata


if __name__ == '__main__':
    #todo interzone not there...
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\number_of_steps"
    if True: # just a way to stop this from running
        summary = run_scenarios(outdir)
    else:
        summary = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\number_of_steps\all_nsteps_summary.csv", index_col=0)

    # generate subset
    keep_cols = [
        'd10_s0_L0_excld_without_pc5pa00_preferred_tar',
        'd10_s5_L0_excld_without_pc5pa00_preferred_tar',
        'd10_s5_L2_excld_without_pc5pa00_preferred_tar',
        'd25_s5_L2_excld_without_pc5pa00_preferred_tar',
        'd25_s5_L2_incld_without_pc5pa00_preferred_tar',
        'd10_s0_L0_incld_with_pc5pa00_preferred_tar',
        'd10_s5_L0_excld_with_pc5pa00_preferred_tar',
        'd10_s0_L0_incld_without_pc5pa00_alternate_tar',
        'd10_s5_L2_incld_without_pc5pa00_alternate_tar',
        'd10_s5_L2_incld_with_pc5pa00_alternate_tar',
    ]
    cols = deepcopy(summary.columns)
    summary.loc['dairy_red',:] = [e.split('_')[0].replace('d','') for e in cols]
    summary.loc['sbda_red',:] = [e.split('_')[1].replace('s','') for e in cols]
    summary.loc['lifestyle_red',:] = [e.split('_')[2].replace('L','') for e in cols]
    summary.loc['pa',:] = [e.split('_')[4].replace('without','pc5pa').replace('with','Pc5pa00') for e in cols]
    summary.loc['low drain',:] = [e.split('_')[3] for e in cols]
    summary.loc['target',:] = [e.split('_')[-2] for e in cols]
    newrows = ['dairy_red','sbda_red','lifestyle_red','pa','target', 'low drain']
    summary.loc[newrows+order,keep_cols].to_csv(os.path.join(outdir, 'a_select_scenarios.csv'))
    print('done')
