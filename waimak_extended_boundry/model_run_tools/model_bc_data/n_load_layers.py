# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/04/2018 11:20 AM
"""

from __future__ import division
import env
from waimak_extended_boundry import smt
import numpy as np
import os
import geopandas as gpd
import shutil
import netCDF4 as nc

def get_gmp_con_layer(recalc=False): #todo test
    """
    'recharge concentration under good managment practices GMP'
    :param recalc:
    :return:
    """
    data = nc.Dataset(os.path.join(env.sdp_required,'N_layers.nc'))
    if (not recalc):
        outdata = np.array(data.variables['gmp_n_conc'])
        return outdata

    raise NotImplementedError ('below is left for documentation only')
    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_new_cmp_con(recalc=False): #todo test
    """
    'recharge concentration under current managment practices'
    :param recalc:
    :return:
    """
    # should be similar to the get original cmp layer, but I am unsure how the re-sampling was handled for that layer
    data = nc.Dataset(os.path.join(env.sdp_required,'N_layers.nc'))
    if (not recalc):
        outdata = np.array(data.variables['cmp_con'])
        return outdata

    raise NotImplementedError('below is left for documentation only')
    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_cmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_pc5pa_additonal_load(recalc=False): # todo test
    """
    get the additional load from pc5pa rules
    :param recalc:
    :return:
    """
    data = nc.Dataset(os.path.join(env.sdp_required,'N_layers.nc'))
    if (not recalc):
        outdata = np.array(data.variables['pc5pa_additional_n_load'])
        return outdata


    raise NotImplementedError('below is left for documentation only')
    n_load_path = r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb"
    layer = 'BAU_PC5_diff_woIrrig_180315'

    outdata = smt.geodb_to_model_array(n_load_path, layer, attribute='BAU_PC5_diff', alltouched=True,
                                       area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_pc5pa_additonal_con(recalc=False):# todo test
    """
    fraction of the additional load with new pc5PA rules
    :param recalc:
    :return:
    """
    data = nc.Dataset(os.path.join(env.sdp_required,'N_layers.nc'))
    if (not recalc):
        outdata = np.array(data.variables['pc5pa_additional_n_con_per'])
        return outdata

    raise NotImplementedError('below is left for documentation only')
    n_load_path = r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb"
    layer = 'BAU_PC5_diff_woIrrig_180315'

    outdata = smt.geodb_to_model_array(n_load_path, layer, attribute='pc5pa_frac_increase', alltouched=True,
                                       area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 1
    np.savetxt(pickle_path, outdata)
    return outdata


def get_gmp_load_raster(recalc=False): #todo test
    """
    array of N load for good management practices (GMP)
    :param recalc:
    :return:
    """
    data = nc.Dataset(os.path.join(env.sdp_required,'N_layers.nc'))
    if (not recalc):
        outdata = np.array(data.variables['gmp_load'])
        return outdata

    raise NotImplementedError('below is left for documentation only')
    n_load_path = os.path.join(env.sdp_required, "NloadLayers/CMP_GMP_PointSources290118.shp")
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nload_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_gmp_plus_con_layer_by_landuse(
        add_half_pc5pa=False,
        exclude_ashley=False,
        Arable=None,
        DairyFarm=None,
        DairySupport=None,
        ForestTussock=None,
        Horticulture=None,
        Lifestyle=None,
        NotFarm=None,
        OtherinclGolf=None,
        Pigs=None,
        SheepBeefDeer=None,
        SheepBeefHill=None,
        Unknow=None):
    """

    get the concentration layer for beyond gmp reductions for a given landuse (defined here and in the shapefile,
    :param exclude_ashley: boolean, if True exlcude the ashley zone area for reductions (see internal shapefile for zone)
    :param add_half_pc5pa: bool if True, add half (to assume reasonable uptake) of teh additional PA N
    # below are land types, reductions are expected in percentage e.g. 20 means 20% reduction from GMP
    :param Arable: all parameters below here are either None or an integer from 0 to 100 which denotes the percentage
                   reduction beyond GMP for the given land use type
    :param DairyFarm:
    :param DairySupport:
    :param ForestTussock:
    :param Horticulture:
    :param Lifestyle:
    :param NotFarm:
    :param OtherinclGolf:
    :param Pigs:
    :param SheepBeefDeer:
    :param SheepBeefHill:
    :param Unknow:
    :return:
    """
    # this takes about 1 minute to run... I'm not pickling because there are too many options, but this
    # could be re-assessed if need to load it lots... better to calculate it once and then make copies...
    base_nload_path = os.path.join(env.sdp_required, "NloadLayers/CMP_GMP_PointSources290118.shp")

    land_types = ('Arable',
                  'DairyFarm',
                  'DairySupport',
                  'ForestTussock',
                  'Horticulture',
                  'Lifestyle',
                  'NotFarm',
                  'OtherinclGolf',
                  'Pigs',
                  'SheepBeefDeer',
                  'SheepBeefHill',
                  'Unknow',)
    base_data = gpd.read_file(base_nload_path)
    base_data.replace({'luscen_cat': {
        'Forest-Tussock': 'ForestTussock',
        'Other-inclGolf': 'OtherinclGolf',
        'Sheep-Beef-Deer': 'SheepBeefDeer',
        'SheepBeef-Hill': 'SheepBeefHill',
    }})
    base_data.loc[:, 'use_n'] = base_data.loc[:, 'nconc_gmp']
    if exclude_ashley:  # for now I'll not worry about it
        raise NotImplementedError
    else:
        for key in land_types:
            if eval(key) is not None:
                base_data.loc[base_data.luscen_cat == key, 'use_n'] *= (1. - eval(key) / 100.)
    temp_shp_path = os.path.join(smt.temp_file_dir, 'temp_ncon_shp', 'temp_ncon_shp.shp')
    if not os.path.exists(os.path.dirname(temp_shp_path)):
        os.makedirs(os.path.dirname(temp_shp_path))
    base_data.to_file(temp_shp_path)

    outdata = smt.shape_file_to_model_array(temp_shp_path, attribute='use_n', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0

    if add_half_pc5pa:
        pa = get_pc5pa_additonal_con()
        pa = (outdata * pa - outdata) / 2
        outdata += pa
    shutil.rmtree(os.path.dirname(temp_shp_path))

    return outdata

def get_gmp_plus_con_layer_by_load_limit(
        n_load_limit,
        n_reduction,
        add_half_pc5pa=False,
        exclude_ashley=False,):
    """
    get the concentration layer for reductions beyond gmp applied to landuse above a n load limit

    :param n_load_limit: array the lower inclusive limit which determines whether or not a parcel will be included
                         in the reductions if nload_limt is longer than 1 the reduction will apply to [lower, next lower)
                         set
    :param n_reduction: percentage (e.g 5 = 5%) array of the amount of reduction to apply limit and reduction are
                        applied positionally
    :param add_half_pc5pa: bool if True, add half (to assume reasonable uptake) of teh additional PA N
    :param exclude_ashley: boolean, if True exlcude the ashley zone area for reductions (see internal shapefile for zone)
    :return:
    """
    # this takes about 1 minute to run... I'm not pickling because there are too many options, but this
    # could be re-assessed if need to load it lots... better to calculate it once and then make copies...
    n_load_limit = np.atleast_1d(n_load_limit)
    n_reduction = np.atleast_1d(n_reduction)
    idx = np.argsort(n_load_limit)
    n_load_limit = n_load_limit[idx]
    n_reduction = n_reduction[idx]

    assert n_load_limit.shape == n_reduction.shape, 'n_load_limit and n_reduction must have the same shape'

    base_nload_path = os.path.join(env.sdp_required, "NloadLayers/CMP_GMP_PointSources290118.shp")

    base_data = gpd.read_file(base_nload_path)
    base_data.loc[:, 'use_n'] = base_data.loc[:, 'nconc_gmp']
    if exclude_ashley:  # for now I'll not worry about it
        raise NotImplementedError
    else:
        for i in range(len(n_load_limit)):
            lower = n_load_limit[i]
            if i != len(n_load_limit) - 1: # not the last entry
                upper = n_load_limit[i+1]
            else:
                upper = 999999999999999999999999 # everything is included in the upper
            base_data.loc[(base_data.nload_gmp >= lower) & (base_data.nload_gmp < upper), 'use_n'] *= (1. - n_reduction[i] / 100.)
    temp_shp_path = os.path.join(smt.temp_file_dir, 'temp_ncon_shp', 'temp_ncon_shp.shp')
    if not os.path.exists(os.path.dirname(temp_shp_path)):
        os.makedirs(os.path.dirname(temp_shp_path))
    base_data.to_file(temp_shp_path)

    outdata = smt.shape_file_to_model_array(temp_shp_path, attribute='use_n', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0

    if add_half_pc5pa:
        pa = get_pc5pa_additonal_con()
        pa = (outdata * pa - outdata) / 2
        outdata += pa
    shutil.rmtree(os.path.dirname(temp_shp_path))

    return outdata


if __name__ == '__main__':
    org = get_gmp_con_layer()
    new = get_gmp_plus_con_layer_by_load_limit(n_load_limit = [10,15,25],
                                                n_reduction = [5,10,20])
    print ('done')