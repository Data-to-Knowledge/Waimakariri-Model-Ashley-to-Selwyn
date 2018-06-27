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
from percentage_reduction_maps import gen_stream_targets, gen_well_targets, get_current_pathway_n


use_soil_mapper = {
            'dairy_lowdrain':1,
            'sbda_lowdrain':1,
            'dairy_highdrain':0,
            'sbda_highdrain':0,
}
use_luse_mapper = {
            'dairy_lowdrain':1,
            'sbda_lowdrain':2,
            'dairy_highdrain':1,
            'sbda_highdrain':2,
}

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
            idx = paths.str.contains(site.replace('wdc_','')) & paths.str.contains('wdc_use_mix')
        else:
            idx = paths.str.contains(site) & ~paths.str.contains('wdc_use_mix')

        outdata[site] = paths.loc[idx].iloc[0]

    return outdata


def get_average_n(recalc=False):
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
                idx = (temp_n.use_ltype==use_luse_mapper[lclass]) & (temp_n.use_soil == use_soil_mapper[lclass])
            total_area = temp_n.loc[idx, 'use_area'].sum()
            if total_area ==0:
                n = 0
            else:
                n = (temp_n.loc[idx, 'use_area'] * temp_n.loc[idx, 'nconc_gmp']).sum() / total_area
            temp_out[lclass] = n
        outdata[site] = temp_out
    return outdata


def get_fractions():
    # todo get the fraction for each component for each site and landuse
    # missing fraction denoted as other
    raise NotImplementedError


def number_of_steps(step_reductions, pa_00, target_scheme):
    """
    calculate the number of steps to reach a selected target scheme
    :param step_reductions: a dictionary of the percentage reduction to apply at each step
                            keys: 'dairy_lowdrain',
                            'sbda_lowdrain',
                            'dairy_highdrain',
                            'sbda_highdrain',
                            'lifestyle'
    :param pa_00: boolean if True, set teh pa rules to 0,0
    :param target_scheme: one of: ['preferred', 'alternate']
    :return:
    """
    assert isinstance(step_reductions, dict)
    assert set(step_reductions.keys()) == lclasses

    targets = gen_well_targets(target_scheme)
    targets.update(gen_stream_targets(target_scheme))
    mode = {k: '50%' for k in targets.keys()}
    current_paths = get_current_pathway_n(mode=mode, conservative_zones='use_mix', from_mt3d_runs=True,
                                          mt3d_add_pa=not pa_00)

    average_n = get_average_n()
    fractions = get_fractions()
    outdata = pd.DataFrame(index=targets.keys())
    for site in targets.keys():
        # calculate the average n concentration for the missing components (check this)
        n_proj = current_paths.loc[site]
        n_proj_temp = current_paths.loc[site]
        for lclass in lclasses:
            n_proj_temp += - average_n[site][lclass] * fractions[site][lclass]

        other_n_con = n_proj_temp/fractions[site]['other']
        outdata.loc[site,'other_n_conc'] = other_n_con

        # apply reduction, calculate new concentration
        new_n = other_n_con * fractions[site]['other']
        for lclass in lclasses:
            new_n += average_n[site][lclass] * step_reductions[lclass] * fractions[site][lclass]
        outdata.loc[site, 'step_1_con'] = new_n

        # from new concentration calculate the number of steps to the target
        dif = n_proj - new_n
        num_steps = (n_proj - targets[site])/dif
        outdata.loc[site, 'nsteps'] = num_steps

    return outdata #todo check this

if __name__ == '__main__':
    test = get_average_n()
    print('done')
