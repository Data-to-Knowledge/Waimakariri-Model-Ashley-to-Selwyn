# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/04/2018 11:20 AM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import numpy as np
import os
import flopy
import geopandas as gpd
import shutil


def get_gmp_con_layer(recalc=False):
    pickle_path = "{}/gmp_n_conc.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_new_cmp_con(recalc=False):
    # should be similar to the get original cmp layer, but I am unsure how the re-sampling was handled for that layer
    pickle_path = "{}/cmp_n_conc_resampled.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_cmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_pc5pa_additonal_load(recalc=False):
    """
    get the additional load from pc5pa
    :param recalc:
    :return:
    """
    pickle_path = "{}/pc5pa_additional_n_load.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb"
    layer = 'BAU_PC5_diff_woIrrig_180315'

    outdata = smt.geodb_to_model_array(n_load_path, layer, attribute='BAU_PC5_diff', alltouched=True,
                                       area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_pc5pa_additonal_con(recalc=False):
    """
    get the additional concentration increase as fraction from pc5pa
    :param recalc:
    :return:
    """
    pickle_path = "{}/pc5pa_additional_n_con_per.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb"
    layer = 'BAU_PC5_diff_woIrrig_180315'

    outdata = smt.geodb_to_model_array(n_load_path, layer, attribute='pc5pa_frac_increase', alltouched=True,
                                       area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 1
    np.savetxt(pickle_path, outdata)
    return outdata


def get_gmp_load_raster(recalc=False):
    pickle_path = "{}/gmp_n_load.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nload_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_orginal_cmp_layer():
    # the layer that all of brioch's runs were done with
    rch_path = env.gw_met_data(r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_"
                               r"Ash0\nconc_cmp_200m.ref")

    return flopy.utils.Util2d.load_txt((smt.rows, smt.cols), rch_path, float, '(FREE)')


def get_gmp_plus_con_layer(
        add_half_pc5pa=False,
        exclude_ashley=True,
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

    :param exclude_ashley: boolean, if True exlcude the ashley zone area for reductions (see internal shapefile for zone)
    :param add_half_pc5pa: bool if True, add half (to assume reasonable uptake) of teh additional PA N
    # below are land types, reductions are expected in percentage e.g. 20 means 20% reduction from GMP
    :param Arable:
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
    base_nload_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                              'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')

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


if __name__ == '__main__':
    org = get_gmp_con_layer()
    test = get_gmp_plus_con_layer(False,False,DairyFarm=25,DairySupport=25)
    if True:
        gmp = get_gmp_con_layer()
        cmp = get_new_cmp_con()
        temp = gmp / cmp
        smt.array_to_raster(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\interzone_n_results\gmpcon_over_cmpcon.tif",
            temp, no_flow_layer=0)
