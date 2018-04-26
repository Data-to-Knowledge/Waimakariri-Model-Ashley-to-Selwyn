# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 23/04/2018 9:02 AM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
import netCDF4 as nc
import itertools
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.source_delmination.interzone_source_delineation import \
    base_receptors_path
import os
import psutil
import geopandas as gpd

p = psutil.Process(os.getpid())
# set to lowest priority, this is windows only, on Unix use ps.nice(19)
p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)

layer_names = {
    1: 'Riccaton',
    3: 'Linwood',
    5: 'Burwood',
    7: 'Wainoni',
    8: 'deep unnamed',
    9: 'deep unnamed',
    10: 'who knows whats happening here'

}

scenario_paths = {
    'cmp': r"C:\mh_waimak_model_data\mednload_ucn.nc",
    # this one is compressed! env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\mednload_unc.nc"),
    'gmp': env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn.nc"),
    'interzone_8kgha': env.gw_met_data(
        r"mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_8kg_ha_interzone.nc"),
    'chch_8kgha': env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_8kg_ha_chch.nc"),
    'gmp_eyre_mar': env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_eyre_mar_ucn.nc"),
    'interzone_50_red': env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\GMP_mednload_ucn_50_reduc_interzone.nc")
}


def make_shapefiles(outdir):
    # make shapefiles for viewing
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    base_pieces = gpd.read_file(base_receptors_path)
    for site, zones in zone_nums.items():
        temp = base_pieces.loc[np.in1d(base_pieces.zone, zones)]
        temp.loc[:, 'grouper'] = 1
        temp = temp.dissolve(by='grouper')
        temp.to_file(os.path.join(outdir, site + '.shp'))


zone_nums = {
    'north_western_tla': [1, 2],
    'western_chch': [6, 7, 8],
    'north_east': [11, 14, 15, 4],
    'central_chch': [10, 13],
    'southern_chch': [9, 12],
    'full_city': [8, 11, 14, 7, 10, 13, 6, 9, 12]
}


def get_chch_area_zones():
    base_zone = smt.shape_file_to_model_array(base_receptors_path, 'zone', True)
    out_zones = {}
    for zone in zone_nums:
        temp = smt.get_empty_model_grid().astype(bool)
        for num in zone_nums[zone]:
            temp[np.isclose(base_zone, num)] = True
        out_zones[zone] = temp
    return out_zones


def get_interzone_n(scenarios, outpath):
    all_indexes = get_chch_area_zones()
    layers = [[1], [3], [5], [7], [8, 9]]  # lump the big deep layers
    index_names = all_indexes.keys()
    # initalize the outdata multi-index dataframe
    percentiles = (0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)
    pers_names = _np_describe(smt.get_empty_model_grid(), percentiles).index.values
    outdata = pd.DataFrame(index=pd.MultiIndex.from_product(([layer_names[e[0]] for e in layers], index_names),
                                                            names=['aquifer', 'zone']),
                           columns=pd.MultiIndex.from_product((scenarios, pers_names),
                                                              names=['scenario', 'stat']), dtype=float)

    for lay, idx_name, scen in itertools.product(layers, index_names, scenarios):
        print(lay, idx_name, scen)
        layer_name = layer_names[lay[0]]
        raw_data = get_raw_model_results(scen, lay, all_indexes[idx_name])
        stocastic_data = apply_nload_uncertainty(raw_data, scen)
        t = _np_describe(np.random.choice(stocastic_data, 100000), percentiles)
        for per in pers_names:
            outdata.loc[(layer_name, idx_name), (scen, per)] = t[per]
    outdata.to_csv(outpath)


def get_raw_model_results(scenario, layers, index):
    nsmc_nums = get_stocastic_set(False)
    data = nc.Dataset(scenario_paths[scenario])
    nsmc_idx = np.in1d(data.variables['nsmc_num'][:], nsmc_nums)
    out_data = np.array(data.variables['mednload'][nsmc_idx, layers])[:, :, index].flatten()

    return out_data


def apply_nload_uncertainty(raw_data, scenario):
    if 'gmp' in scenario:
        modifiers = np.loadtxt(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and "
                               r"results\ex_bd_va\n_results\interzone\stocastic_n_interzone\without_trans\nload_cmp_"
                               r"interzone\raw_data\conservative_interzone.txt")[
                    0:1000]  # cmp because it's referenced to cmp
    elif 'cmp' in scenario:
        modifiers = np.loadtxt(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and "
                               r"results\ex_bd_va\n_results\interzone\stocastic_n_interzone\with_trans\nconc_cmp_"
                               r"interzone\raw_data\conservative_interzone.txt")[0:1000]
    else:
        return raw_data
    max_size = 2000 * 165 * 2
    raw_data = raw_data.astype(np.float32)
    if len(raw_data) > max_size:
        raw_data = np.random.choice(raw_data, (max_size))  # to get around run size errors
    all_n = (raw_data[:, np.newaxis] * modifiers[np.newaxis, :].astype(np.float32)).flatten()
    return all_n


def _np_describe(ndarray, percentiles=(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)):
    """
    pull out the percentiles and put it in a pd.Series
    :param ndarray: a 1 or more dimentional numpy array
    :param percentiles: percentiles to return
    :return: pd.Series
    """
    outdata = pd.Series()
    outdata.loc['mean'] = np.nanmean(ndarray)
    outdata.loc['std'] = np.nanstd(ndarray)
    outdata.loc['min'] = np.nanmin(ndarray)
    for i in percentiles:
        outdata.loc['{}%'.format(int(i * 100))] = np.nanpercentile(ndarray, i * 100)
    outdata.loc['max'] = np.nanmax(ndarray)
    return outdata


def make_shpfile_with_data(n_data_path, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    names = {
        u'chch_8kgha_50%': 'chch50',
        u'chch_8kgha_95%': 'chch95',
        u'interzone_8kgha_50%': 'inter50',
        u'interzone_8kgha_95%': 'inder95',
        'interzone_50_red_50%': 'in50_50',
        'interzone_50_red_95%': 'in50_95',
        u'gmp_50%': 'gmp_50',
        u'gmp_95%': 'gmp_95',
        u'gmp_eyre_mar_50%': 'mar50',
        u'gmp_eyre_mar_95%': 'mar95',
        u'cmp_50%': 'cmp_50',
        u'cmp_95%': 'cmp_95',
        u'chch_8kgha_plus_interzone_8kgha_50%': '8_8_50',
        u'chch_8kgha_plus_interzone_8kgha_95%': '8_8_95',
        u'chch_8kgha_plus_gmp_50%': '8_gmp50',
        u'chch_8kgha_plus_gmp_95%': '8_gmp95',
        u'chch_8kgha_plus_gmp_eyre_mar_50%': '8_mar50',
        u'chch_8kgha_plus_gmp_eyre_mar_95%': '8_mar95',
        u'chch_8kgha_plus_cmp_50%': '8_cmp50',
        u'chch_8kgha_plus_cmp_95%': '8_cmp95',
        'chch_8kgha_plus_interzone_50_red_50%': '8_in50_50',
        'chch_8kgha_plus_interzone_50_red_95%': '8_in50_95'}

    layers = [[1], [3], [5], [7], [8, 9]]
    all_data = pd.read_csv(
        n_data_path,
        index_col=[0, 1], header=[0, 1])
    for lay in layers:
        lay_name = layer_names[lay[0]]

        base_pieces = gpd.read_file(base_receptors_path)
        for site, zones in zone_nums.items():
            temp = base_pieces.loc[np.in1d(base_pieces.zone, zones)]
            temp.loc[:, 'grouper'] = 1
            temp = temp.dissolve(by='grouper')

            for scen, per in itertools.product(scenario_paths.keys(), ['50%', '95%']):
                temp.loc[1, '{}_{}'.format(scen, per)] = all_data.loc[(lay_name, site), (scen, per)]

            for scen in ['interzone_8kgha', 'gmp', 'gmp_eyre_mar', 'cmp', 'interzone_50_red']:
                for per in ['50%', '95%']:
                    temp.loc[1, 'chch_8kgha_plus_{}_{}'.format(scen, per)] = temp.loc[1, 'chch_8kgha_{}'.format(per)] + \
                                                                             temp.loc[1, '{}_{}'.format(scen, per)]
            temp = temp.rename(columns=names)
            temp.to_file(os.path.join(outdir, '{}_{}.shp'.format(site, lay_name)))


if __name__ == '__main__':
    if False:
        make_shapefiles(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex"
                        r"_bd_va\n_results\interzone_n_results\chch_wm_receptor_shapes")

    if False:
        get_interzone_n(scenario_paths.keys(),
                        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                        r"s\ex_bd_va\n_results\interzone_n_results\n_data_v2.csv")
        data = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                           r"s\ex_bd_va\n_results\interzone_n_results\n_data_v2.csv",
                           index_col=[0,1], header=[0,1])
        data = data.reorder_levels(['stat','scenario'],axis=1)
        data['50%'].to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                           r"s\ex_bd_va\n_results\interzone_n_results\n_data_50ths_v2.csv")
        data['95%'].to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                           r"s\ex_bd_va\n_results\interzone_n_results\n_data_95ths_v2.csv")
        data['mean'].to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and result"
                            r"s\ex_bd_va\n_results\interzone_n_results\n_data_means_v2.csv")
    if True:
        make_shpfile_with_data(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\interzone_n_results\n_data_v2.csv",
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\interzone_n_results\receptors_with_data_inc_50_red")
