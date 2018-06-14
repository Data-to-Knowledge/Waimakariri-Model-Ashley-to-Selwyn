"""
Author: matth
Date Created: 20/06/2017 11:57 AM
"""

from __future__ import division
from core import env
import flopy_mh as flopy
from users.MH.Waimak_modeling.model_tools.well_values import get_race_data, get_nwai_wells
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.supporting_data_path import sdp
from core.ecan_io import rd_sql, sql_db
import os
import pickle
import geopandas as gpd
from copy import deepcopy
from warnings import warn


def create_wel_package(m, wel_version):
    """
    create and add the well package
    :param m: a flopy model instance
    :param wel_version: which version of wells to use
    :return:
    """
    wel = flopy.modflow.mfwel.ModflowWel(m,
                                         ipakcb=740,
                                         stress_period_data={
                                             0: smt.convert_well_data_to_stresspd(get_wel_spd(wel_version))},
                                         options=['AUX IFACE'],  # next time don't include this unless I use it
                                         unitnumber=709)


def get_wel_spd(version, recalc=False):
    """
    get the well data
    :param version: which well version to use
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return: pd.dataframe
    """
    if version == 1:
        outdata = _get_wel_spd_v1(recalc, sub_version=1)
    elif version == 0:
        outdata = _get_wel_spd_v1(recalc, sub_version=0)
    elif version == 3:
        outdata = _get_wel_spd_v3(recalc)
    else:
        raise ValueError('unexpected version: {}'.format(version))
    return outdata


def _get_wel_spd_v1(recalc=False, sub_version=1):
    """
    version 1, which is deprecated but uses use data and other estimates for the waimakariri zone
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :param sub_version: passed to get swai wells
    :return: pd.DataFrame
    """
    warn('v1 wells are depreciated in the newest itteration of the model')
    pickle_path = '{}/well_spd.p'.format(smt.pickle_dir)
    if os.path.exists(pickle_path) and not recalc and sub_version != 0:
        well_data = pickle.load(open(pickle_path))
        return well_data

    races = get_race_data()

    elv_db = smt.calc_elv_db()
    for site in races.index:
        races.loc[site, 'row'], races.loc[site, 'col'] = smt.convert_coords_to_matix(races.loc[site, 'x'],
                                                                                     races.loc[site, 'y'])
    races['zone'] = 'n_wai'
    races = races.set_index('well')

    n_wai_wells = get_nwai_wells()
    for site in n_wai_wells.index:
        x, y, z = n_wai_wells.loc[site, ['x', 'y', 'z']]
        temp = smt.convert_coords_to_matix(x, y, z, elv_db=elv_db)
        n_wai_wells.loc[site, 'layer'], n_wai_wells.loc[site, 'row'], n_wai_wells.loc[site, 'col'] = temp
    n_wai_wells['zone'] = 'n_wai'
    n_wai_wells['cwms'] = 'waimak'
    n_wai_wells = n_wai_wells.set_index('well')

    s_wai_wells = _get_s_wai_wells(
        sub_version)  # there are some s_wai wells which do not have data in wells, but do in consents file fix if bored
    temp = smt.get_well_postions(np.array(s_wai_wells.index), one_val_per_well=True, raise_exct=False)
    s_wai_wells['layer'], s_wai_wells['row'], s_wai_wells['col'] = temp
    no_flow = smt.get_no_flow()
    for i in s_wai_wells.index:
        layer, row, col = s_wai_wells.loc[i, ['layer', 'row', 'col']]
        if any(pd.isnull([layer, row, col])):
            continue
        if no_flow[layer, row, col] == 0:  # get rid of non-active wells
            s_wai_wells.loc[i, 'layer'] = np.nan

    s_wai_wells = s_wai_wells.dropna(subset=['layer', 'row', 'col'])

    s_wai_rivers = _get_s_wai_rivers().set_index('well')

    all_wells = pd.concat((races, n_wai_wells, s_wai_wells, s_wai_rivers))

    for i in all_wells.index:
        row, col = all_wells.loc[i, ['row', 'col']]
        x, y = smt.convert_matrix_to_coords(row, col)
        all_wells.loc[i, 'mx'] = x
        all_wells.loc[i, 'my'] = y

    # check wells in correct aquifer
    aq_to_layer = {'Avonside Formation': 0,
                   'Springston Formation': 0,
                   'Christchurch Formation': 0,
                   'Riccarton Gravel': 1,
                   'Bromley Formation': 2,
                   'Linwood Gravel': 3,
                   'Heathcote Formation': 4,
                   'Burwood Gravel': 5,
                   'Shirley Formation': 6,
                   'Wainoni Gravel': 7}
    leapfrog_aq = gpd.read_file("{}/m_ex_bd_inputs/shp/layering/gis_aq_name_clipped.shp".format(smt.sdp))
    leapfrog_aq = leapfrog_aq.set_index('well')
    leapfrog_aq.loc[:, 'use_aq_name'] = leapfrog_aq.loc[:, 'aq_name']
    leapfrog_aq.loc[leapfrog_aq.use_aq_name.isnull(), 'use_aq_name'] = leapfrog_aq.loc[
        leapfrog_aq.use_aq_name.isnull(), 'aq_name_gi']

    for well in all_wells.index:
        try:
            all_wells.loc[well, 'aquifer_in_confined'] = aq = leapfrog_aq.loc[well, 'use_aq_name']
            all_wells.loc[well, 'layer_by_aq'] = aq_to_layer[aq]
        except KeyError:
            pass

    all_wells.loc[:, 'layer_by_depth'] = all_wells.loc[:, 'layer']
    all_wells.loc[all_wells.layer_by_aq.notnull(), 'layer'] = all_wells.loc[
        all_wells.layer_by_aq.notnull(), 'layer_by_aq']

    # move wells that fall on other boundry conditions north of waimak (or in constant head)
    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2.shp".format(smt.sdp))
    overlap = overlap.set_index('index')

    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_k']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2part2.shp".format(smt.sdp))
    overlap = overlap.set_index('Field1')
    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_layer']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    # add little rakaia flux which will be parameterized via pest in two groups upper flux is north of SH1, lower is coastal of SH1
    temp = smt.model_where(np.isfinite(
        smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/little_rakaia_boundry_wells.shp".format(smt.sdp),
                                      'Id', True)))
    all_llrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        llrf = pd.DataFrame(index=['llrz_flux{:04d}'.format(e) for e in range(i * len(temp), (i + 1) * len(temp))],
                            columns=all_wells.keys())
        llrf.loc[:, 'row'] = np.array(temp)[:, 0]
        llrf.loc[:, 'col'] = np.array(temp)[:, 1]
        llrf.loc[:, 'layer'] = i
        llrf.loc[:, 'flux'] = -9999999  # identifier flux, parameterised in pest
        llrf.loc[:, 'type'] = 'llr_boundry_flux'
        llrf.loc[:, 'zone'] = 's_wai'
        all_llrf = pd.concat((all_llrf, llrf))

    up_temp = smt.model_where(
        np.isfinite(smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/upper_lRZF.shp".format(smt.sdp),
                                                  'Id', True)))

    all_ulrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        ulrf = pd.DataFrame(
            index=['ulrz_flux{:04d}'.format(e) for e in range(i * len(up_temp), (i + 1) * len(up_temp))],
            columns=all_wells.keys())
        ulrf.loc[:, 'row'] = np.array(up_temp)[:, 0]
        ulrf.loc[:, 'col'] = np.array(up_temp)[:, 1]
        ulrf.loc[:, 'layer'] = i
        ulrf.loc[:, 'flux'] = -8888888  # identifier flux, parameterised in pest
        ulrf.loc[:, 'type'] = 'ulr_boundry_flux'
        ulrf.loc[:, 'zone'] = 's_wai'
        all_ulrf = pd.concat((all_ulrf, ulrf))

    swai_races = get_s_wai_races()
    all_wells = pd.concat((all_wells, swai_races, all_llrf, all_ulrf))

    all_wells = all_wells.loc[~((all_wells.duplicated(subset=['row', 'col', 'layer'], keep=False)) &
                                (all_wells.type.str.contains('lr_boundry_flux')))]
    all_wells = add_use_type(all_wells)  # any well that has irrigation/stockwater in it's uses is considered irrigation
    if sub_version != 0:
        pickle.dump(all_wells, open(pickle_path, 'w'))
    return all_wells


def _get_wel_spd_v3(recalc=False, sub_version=1):
    """
    all wells derived from mikes usage estimates I may pull down some of the WDC WS wells this was used in teh model as
    of 20/10/2017
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :param sub_version: passed to get all wells
    :return: pd.DataFrame
    """
    pickle_path = '{}/well_spd_v3.p'.format(smt.pickle_dir)
    if os.path.exists(pickle_path) and not recalc and sub_version != 0:
        well_data = pickle.load(open(pickle_path))
        return well_data

    races = get_race_data()

    elv_db = smt.calc_elv_db()
    for site in races.index:
        races.loc[site, 'row'], races.loc[site, 'col'] = smt.convert_coords_to_matix(races.loc[site, 'x'],
                                                                                     races.loc[site, 'y'])
    races['zone'] = 'n_wai'
    races = races.set_index('well')
    all_wai_wells = _get_all_wai_wells()
    temp = smt.get_well_postions(np.array(all_wai_wells.index), one_val_per_well=True, raise_exct=False)
    all_wai_wells['layer'], all_wai_wells['row'], all_wai_wells['col'] = temp
    no_flow = smt.get_no_flow()
    for i in all_wai_wells.index:
        layer, row, col = all_wai_wells.loc[i, ['layer', 'row', 'col']]
        if any(pd.isnull([layer, row, col])):
            continue
        if no_flow[layer, row, col] == 0:  # get rid of non-active wells
            all_wai_wells.loc[i, 'layer'] = np.nan

    all_wai_wells = all_wai_wells.dropna(subset=['layer', 'row', 'col'])

    s_wai_rivers = _get_s_wai_rivers().set_index('well')

    all_wells = pd.concat((races, all_wai_wells, s_wai_rivers))

    for i in all_wells.index:
        row, col = all_wells.loc[i, ['row', 'col']]
        x, y = smt.convert_matrix_to_coords(row, col)
        all_wells.loc[i, 'mx'] = x
        all_wells.loc[i, 'my'] = y

    # check wells in correct aquifer
    aq_to_layer = {'Avonside Formation': 0,
                   'Springston Formation': 0,
                   'Christchurch Formation': 0,
                   'Riccarton Gravel': 1,
                   'Bromley Formation': 2,
                   'Linwood Gravel': 3,
                   'Heathcote Formation': 4,
                   'Burwood Gravel': 5,
                   'Shirley Formation': 6,
                   'Wainoni Gravel': 7}
    leapfrog_aq = gpd.read_file("{}/m_ex_bd_inputs/shp/layering/gis_aq_name_clipped.shp".format(smt.sdp))
    leapfrog_aq = leapfrog_aq.set_index('well')
    leapfrog_aq.loc[:, 'use_aq_name'] = leapfrog_aq.loc[:, 'aq_name']
    leapfrog_aq.loc[leapfrog_aq.use_aq_name.isnull(), 'use_aq_name'] = leapfrog_aq.loc[
        leapfrog_aq.use_aq_name.isnull(), 'aq_name_gi']

    for well in all_wells.index:
        try:
            all_wells.loc[well, 'aquifer_in_confined'] = aq = leapfrog_aq.loc[well, 'use_aq_name']
            all_wells.loc[well, 'layer_by_aq'] = aq_to_layer[aq]
        except KeyError:
            pass

    all_wells.loc[:, 'layer_by_depth'] = all_wells.loc[:, 'layer']
    all_wells.loc[all_wells.layer_by_aq.notnull(), 'layer'] = all_wells.loc[
        all_wells.layer_by_aq.notnull(), 'layer_by_aq']

    # move wells that fall on other boundry conditions north of waimak (or in constant head)
    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2.shp".format(smt.sdp))
    overlap = overlap.set_index('index')

    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_k']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2part2.shp".format(smt.sdp))
    overlap = overlap.set_index('Field1')
    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_layer']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustmentpart3.shp".format(smt.sdp))
    overlap = overlap.set_index('Field1')
    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_layer']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    # note there are some overlaps remaining, but it's probably not a huge problem most are races

    # add little rakaia flux which will be parameterized via pest in two groups upper flux is north of SH1, lower is coastal of SH1
    temp = smt.model_where(np.isfinite(
        smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/little_rakaia_boundry_wells.shp".format(smt.sdp),
                                      'Id', True)))
    all_llrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        llrf = pd.DataFrame(index=['llrz_flux{:04d}'.format(e) for e in range(i * len(temp), (i + 1) * len(temp))],
                            columns=all_wells.keys())
        llrf.loc[:, 'row'] = np.array(temp)[:, 0]
        llrf.loc[:, 'col'] = np.array(temp)[:, 1]
        llrf.loc[:, 'layer'] = i
        llrf.loc[:, 'flux'] = -9999999  # identifier flux, parameterised in pest
        llrf.loc[:, 'type'] = 'llr_boundry_flux'
        llrf.loc[:, 'zone'] = 's_wai'
        all_llrf = pd.concat((all_llrf, llrf))

    up_temp = smt.model_where(
        np.isfinite(smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/upper_lRZF.shp".format(smt.sdp),
                                                  'Id', True)))

    all_ulrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        ulrf = pd.DataFrame(
            index=['ulrz_flux{:04d}'.format(e) for e in range(i * len(up_temp), (i + 1) * len(up_temp))],
            columns=all_wells.keys())
        ulrf.loc[:, 'row'] = np.array(up_temp)[:, 0]
        ulrf.loc[:, 'col'] = np.array(up_temp)[:, 1]
        ulrf.loc[:, 'layer'] = i
        ulrf.loc[:, 'flux'] = -8888888  # identifier flux, parameterised in pest
        ulrf.loc[:, 'type'] = 'ulr_boundry_flux'
        ulrf.loc[:, 'zone'] = 's_wai'
        all_ulrf = pd.concat((all_ulrf, ulrf))

    swai_races = get_s_wai_races()
    all_wells = pd.concat((all_wells, swai_races, all_llrf, all_ulrf))

    all_wells = all_wells.loc[~((all_wells.duplicated(subset=['row', 'col', 'layer'], keep=False)) &
                                (all_wells.type.str.contains('lr_boundry_flux')))]
    all_wells = add_use_type(all_wells)  # any well that has irrigation/stockwater in it's uses is considered irrigation
    if sub_version != 0:
        pickle.dump(all_wells, open(pickle_path, 'w'))
    return all_wells


def _get_2014_2015_waimak_usage():
    """
    get the 2014-2015 usage for the waimakariri zone
    :return: pd.DataFrame
    """
    mike = pd.read_hdf("{}/m_ex_bd_inputs/sd_est_all_mon_vol.h5".format(smt.sdp))
    mike = mike.loc[(mike.time >= pd.datetime(2014, 7, 1)) & (mike.take_type == 'Take Groundwater')]
    mike.loc[:, 'd_in_m'] = mike.time.dt.daysinmonth
    data = mike.groupby('wap').aggregate(
        {'usage_est': np.sum, 'crc': ','.join, 'd_in_m': np.sum, 'mon_allo_m3': np.sum})
    data.loc[:, 'flux'] = data.loc[:, 'usage_est'] / (mike.time.max() - pd.datetime(2014, 6, 30)).days
    data.loc[:, 'cav_flux'] = data.loc[:, 'mon_allo_m3'] / (mike.time.max() - pd.datetime(2014, 6, 30)).days

    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.set_index('WELL_NO')
    out_data = pd.merge(data, pd.DataFrame(well_details.loc[:, 'WMCRZone']), left_index=True, right_index=True)
    out_data = out_data.loc[np.in1d(out_data.WMCRZone, [4])]
    out_data.loc[:, 'cwms'] = out_data.loc[:, 'WMCRZone'].replace({7: 'chch', 8: 'selwyn', 4: 'waimak'})
    out_data = out_data.drop('WMCRZone', axis=1)

    out_data['type'] = 'well'
    out_data = add_use_type(out_data)
    idx = (out_data.cwms == 'waimak') & (out_data.use_type == 'other')
    out_data.loc[idx, 'flux'] = out_data.loc[
                                    idx, 'cav_flux'] * 0.25  # this comes from average of the WDC CAV vs usage made before my time I also confirmed with colin as WDC that this is about right
    idx = out_data.flux > out_data.cav_flux
    out_data.loc[idx, 'flux'] = out_data.loc[idx, 'cav_flux']

    out_data.loc[:, 'flux'] *= -1
    out_data['consent'] = [tuple(e.split(',')) for e in out_data.loc[:, 'crc']]
    out_data = out_data.drop('crc', axis=1)
    out_data = out_data.dropna()

    out_data['zone'] = 'n_wai'
    out_data.index.names = ['well']

    return out_data


def _get_wel_spd_v2(recalc=False, sub_version=1):
    """
    as version 1 but uses the 2014-2015 usage for the waimakariri zone used for the forward runs
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :param sub_version: passed to get s wai wells
    :return: pd.Dataframe
    """
    warn('v2 pumping is for 2014 to 2015 period')
    pickle_path = '{}/well_spd_v2.p'.format(smt.pickle_dir)
    if os.path.exists(pickle_path) and not recalc and sub_version != 0:
        well_data = pickle.load(open(pickle_path))
        return well_data

    races = get_race_data()

    elv_db = smt.calc_elv_db()
    for site in races.index:
        races.loc[site, 'row'], races.loc[site, 'col'] = smt.convert_coords_to_matix(races.loc[site, 'x'],
                                                                                     races.loc[site, 'y'])
    races['zone'] = 'n_wai'
    races = races.set_index('well')

    n_wai_wells = _get_2014_2015_waimak_usage()

    s_wai_wells = _get_s_wai_wells(
        sub_version)  # there are some s_wai wells which do not have data in wells, but do in consents file fix if bored
    ns_wells = pd.concat((n_wai_wells, s_wai_wells))
    temp = smt.get_well_postions(np.array(ns_wells.index), one_val_per_well=True, raise_exct=False)
    ns_wells['layer'], ns_wells['row'], ns_wells['col'] = temp
    no_flow = smt.get_no_flow()
    for i in ns_wells.index:
        layer, row, col = ns_wells.loc[i, ['layer', 'row', 'col']]
        if any(pd.isnull([layer, row, col])):
            continue
        if no_flow[layer, row, col] == 0:  # get rid of non-active wells
            ns_wells.loc[i, 'layer'] = np.nan

    ns_wells = ns_wells.dropna(subset=['layer', 'row', 'col'])

    s_wai_rivers = _get_s_wai_rivers().set_index('well')

    all_wells = pd.concat((races, ns_wells, s_wai_rivers))

    for i in all_wells.index:
        row, col = all_wells.loc[i, ['row', 'col']]
        x, y = smt.convert_matrix_to_coords(row, col)
        all_wells.loc[i, 'mx'] = x
        all_wells.loc[i, 'my'] = y

    # check wells in correct aquifer
    aq_to_layer = {'Avonside Formation': 0,
                   'Springston Formation': 0,
                   'Christchurch Formation': 0,
                   'Riccarton Gravel': 1,
                   'Bromley Formation': 2,
                   'Linwood Gravel': 3,
                   'Heathcote Formation': 4,
                   'Burwood Gravel': 5,
                   'Shirley Formation': 6,
                   'Wainoni Gravel': 7}
    leapfrog_aq = gpd.read_file("{}/m_ex_bd_inputs/shp/layering/gis_aq_name_clipped.shp".format(smt.sdp))
    leapfrog_aq = leapfrog_aq.set_index('well')
    leapfrog_aq.loc[:, 'use_aq_name'] = leapfrog_aq.loc[:, 'aq_name']
    leapfrog_aq.loc[leapfrog_aq.use_aq_name.isnull(), 'use_aq_name'] = leapfrog_aq.loc[
        leapfrog_aq.use_aq_name.isnull(), 'aq_name_gi']

    for well in all_wells.index:
        try:
            all_wells.loc[well, 'aquifer_in_confined'] = aq = leapfrog_aq.loc[well, 'use_aq_name']
            all_wells.loc[well, 'layer_by_aq'] = aq_to_layer[aq]
        except KeyError:
            pass

    all_wells.loc[:, 'layer_by_depth'] = all_wells.loc[:, 'layer']
    all_wells.loc[all_wells.layer_by_aq.notnull(), 'layer'] = all_wells.loc[
        all_wells.layer_by_aq.notnull(), 'layer_by_aq']

    # move wells that fall on other boundry conditions north of waimak (or in constant head)
    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2.shp".format(smt.sdp))
    overlap = overlap.set_index('index')

    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_k']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    overlap = gpd.read_file("{}/m_ex_bd_inputs/shp/overlap_adjustment2part2.shp".format(smt.sdp))
    overlap = overlap.set_index('Field1')
    for well in all_wells.index:
        if not well in overlap.index:
            continue
        all_wells.loc[well, 'layer'] += overlap.loc[well, 'add_layer']
        all_wells.loc[well, 'row'] += overlap.loc[well, 'add_row']
        all_wells.loc[well, 'col'] += overlap.loc[well, 'add_col']

    # add little rakaia flux which will be parameterized via pest in two groups upper flux is north of SH1, lower is coastal of SH1
    temp = smt.model_where(np.isfinite(
        smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/little_rakaia_boundry_wells.shp".format(smt.sdp),
                                      'Id', True)))
    all_llrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        llrf = pd.DataFrame(index=['llrz_flux{:04d}'.format(e) for e in range(i * len(temp), (i + 1) * len(temp))],
                            columns=all_wells.keys())
        llrf.loc[:, 'row'] = np.array(temp)[:, 0]
        llrf.loc[:, 'col'] = np.array(temp)[:, 1]
        llrf.loc[:, 'layer'] = i
        llrf.loc[:, 'flux'] = -9999999  # identifier flux, parameterised in pest
        llrf.loc[:, 'type'] = 'llr_boundry_flux'
        llrf.loc[:, 'zone'] = 's_wai'
        all_llrf = pd.concat((all_llrf, llrf))

    up_temp = smt.model_where(
        np.isfinite(smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/upper_lRZF.shp".format(smt.sdp),
                                                  'Id', True)))

    all_ulrf = pd.DataFrame(columns=all_wells.keys())
    for i in range(smt.layers):
        ulrf = pd.DataFrame(
            index=['ulrz_flux{:04d}'.format(e) for e in range(i * len(up_temp), (i + 1) * len(up_temp))],
            columns=all_wells.keys())
        ulrf.loc[:, 'row'] = np.array(up_temp)[:, 0]
        ulrf.loc[:, 'col'] = np.array(up_temp)[:, 1]
        ulrf.loc[:, 'layer'] = i
        ulrf.loc[:, 'flux'] = -8888888  # identifier flux, parameterised in pest
        ulrf.loc[:, 'type'] = 'ulr_boundry_flux'
        ulrf.loc[:, 'zone'] = 's_wai'
        all_ulrf = pd.concat((all_ulrf, ulrf))

    swai_races = get_s_wai_races()
    all_wells = pd.concat((all_wells, swai_races, all_llrf, all_ulrf))

    all_wells = all_wells.loc[~((all_wells.duplicated(subset=['row', 'col', 'layer'], keep=False)) &
                                (all_wells.type.str.contains('lr_boundry_flux')))]
    all_wells = add_use_type(all_wells)  # any well that has irrigation/stockwater in it's uses is considered irrigation
    if sub_version != 0:
        pickle.dump(all_wells, open(pickle_path, 'w'))
    return all_wells


def _get_s_wai_wells(subversion=1):
    """
    get wells south of the river
    :param subversion: if 0 use mike's allo data (depreciated, but held for histories sake)
                       if 1 use mike's usage estimate
    :return: pd.DataFrame
    """
    if subversion == 1:
        mike = pd.read_hdf("{}/m_ex_bd_inputs/sd_est_all_mon_vol.h5".format(smt.sdp))
        mike = mike.loc[(mike.time >= pd.datetime(2008, 1, 1)) & (mike.take_type == 'Take Groundwater')]
        mike.loc[:, 'd_in_m'] = mike.time.dt.daysinmonth
        data = mike.groupby('wap').aggregate({'usage_est': np.sum, 'crc': ','.join, 'd_in_m': np.sum})
        data.loc[:, 'flux'] = data.loc[:, 'usage_est'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days

        well_details = rd_sql(**sql_db.wells_db.well_details)
        well_details = well_details.set_index('WELL_NO')
        out_data = pd.merge(data, pd.DataFrame(well_details.loc[:, 'WMCRZone']), left_index=True, right_index=True)
        out_data = out_data.loc[np.in1d(out_data.WMCRZone, [7, 8])]
        out_data.loc[:, 'cwms'] = out_data.loc[:, 'WMCRZone'].replace({7: 'chch', 8: 'selwyn'})
        out_data = out_data.drop('WMCRZone', axis=1)
        out_data.loc[:, 'flux'] *= -1


    elif subversion == 0:

        allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp), index_col='crc')

        # option 2
        end_time = pd.datetime(2016, 12, 31)
        start_time = pd.datetime(2008, 1, 1)

        allo2 = allo.loc[np.in1d(allo['cwms'], ['Selwyn - Waihora', 'Christchurch - West Melton']) &
                         (allo['take_type'] == 'Take Groundwater') &
                         (
                         (allo.status_details.str.contains('Terminated')) | allo.status_details.str.contains('Issued'))]

        allo2.loc[:, 'to_date'] = pd.to_datetime(allo2.loc[:, 'to_date'], format='%d/%m/%Y', errors='coerce')
        allo2.loc[:, 'from_date'] = pd.to_datetime(allo2.loc[:, 'from_date'], format='%d/%m/%Y', errors='coerce')
        allo2.loc[allo2.loc[:, 'to_date'] > end_time, 'to_date'] = end_time
        allo2.loc[allo2.loc[:, 'to_date'] < start_time, 'to_date'] = None
        allo2.loc[allo2.loc[:, 'from_date'] < start_time, 'from_date'] = start_time
        allo2.loc[allo2.loc[:, 'from_date'] > end_time, 'from_date'] = None
        idx = (pd.notnull(allo2.loc[:, 'to_date']) & pd.notnull(allo2.loc[:, 'from_date']))
        allo2.loc[idx, 'temp_days'] = (allo2.loc[idx, 'to_date'] - allo2.loc[idx, 'from_date'])
        allo2.loc[:, 'days'] = [e.days for e in allo2.loc[:, 'temp_days']]

        # the below appear to be an interal consents marker... and should not be included here as a replacement consent
        # is active at the same time at the consent with negitive number of days
        allo2.loc[allo2.loc[:, 'days'] < 0, 'days'] = 0

        allo2.loc[:, 'flux'] = allo2.loc[:, 'cav'] / 365 * allo2.loc[:, 'days'] / (end_time - start_time).days * -1
        allo2.loc[allo2.use_type == 'irrigation', 'flux'] *= 0.5

        out_data = allo2.reset_index().groupby('wap').aggregate({'flux': np.sum, 'crc': ','.join})
        out_data.loc[:, 'flux'] *= 0.50  # start with a 50% scaling factor from CAV come back if time
    else:
        raise ValueError('unexpected sub-version')
    out_data['consent'] = [tuple(e.split(',')) for e in out_data.loc[:, 'crc']]
    out_data = out_data.drop('crc', axis=1)
    out_data = out_data.dropna()

    out_data['type'] = 'well'
    out_data['zone'] = 's_wai'
    out_data.index.names = ['well']

    return out_data


def _get_all_wai_wells():
    """
    get's well data from mike's usage for V3
    :return: pd.DataFrame
    """
    # there are some wells where the flux is greater than the CAV;
    # however these are rather minor, moslty in Selwyn, and most could be true.
    mike = pd.read_hdf("{}/m_ex_bd_inputs/sd_est_all_mon_vol.h5".format(smt.sdp))
    mike = mike.loc[(mike.time >= pd.datetime(2008, 1, 1)) & (mike.take_type == 'Take Groundwater')]
    mike.loc[:, 'd_in_m'] = mike.time.dt.daysinmonth
    data = mike.groupby('wap').aggregate(
        {'usage_est': np.sum, 'crc': ','.join, 'd_in_m': np.sum, 'mon_allo_m3': np.sum})
    data.loc[:, 'flux'] = data.loc[:, 'usage_est'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days
    data.loc[:, 'cav_flux'] = data.loc[:, 'mon_allo_m3'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days

    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.set_index('WELL_NO')
    out_data = pd.merge(data, pd.DataFrame(well_details.loc[:, 'WMCRZone']), left_index=True, right_index=True)
    out_data = out_data.loc[np.in1d(out_data.WMCRZone, [4, 7, 8])]
    out_data.loc[:, 'cwms'] = out_data.loc[:, 'WMCRZone'].replace({7: 'chch', 8: 'selwyn', 4: 'waimak'})
    out_data = out_data.drop('WMCRZone', axis=1)

    out_data['type'] = 'well'
    out_data = add_use_type(out_data)

    # set WDC (waimak and other usage) wells to 10% of CAV
    idx = (out_data.cwms == 'waimak') & (out_data.use_type == 'other')
    out_data.loc[idx, 'flux'] = out_data.loc[
                                    idx, 'cav_flux'] * 0.25  # this comes from average of the WDC CAV vs usage made before my time I also confirmed with colin as WDC that this is about right

    out_data.loc[:, 'flux'] *= -1

    out_data['consent'] = [tuple(e.split(',')) for e in out_data.loc[:, 'crc']]
    out_data = out_data.drop('crc', axis=1)
    out_data = out_data.dropna()

    out_data.loc[out_data.cwms == 'waimak', 'zone'] = 'n_wai'
    out_data.loc[~(out_data.cwms == 'waimak'), 'zone'] = 's_wai'
    out_data.index.names = ['well']

    return out_data


def _check_chch_wells():
    """
    some well checks
    :return:
    """
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp), index_col='crc')

    # option 2
    end_time = pd.datetime(2016, 12, 31)
    start_time = pd.datetime(2008, 1, 1)

    allo2 = allo.loc[np.in1d(allo['cwms'], ['Christchurch - West Melton']) &
                     (allo['take_type'] == 'Take Groundwater') &
                     ((allo.status_details.str.contains('Terminated')) | allo.status_details.str.contains('Issued'))]

    allo2.loc[:, 'to_date'] = pd.to_datetime(allo2.loc[:, 'to_date'], format='%d/%m/%Y', errors='coerce')
    allo2.loc[:, 'from_date'] = pd.to_datetime(allo2.loc[:, 'from_date'], format='%d/%m/%Y', errors='coerce')
    allo2.loc[allo2.loc[:, 'to_date'] > end_time, 'to_date'] = end_time
    allo2.loc[allo2.loc[:, 'to_date'] < start_time, 'to_date'] = None
    allo2.loc[allo2.loc[:, 'from_date'] < start_time, 'from_date'] = start_time
    allo2.loc[allo2.loc[:, 'from_date'] > end_time, 'from_date'] = None
    idx = (pd.notnull(allo2.loc[:, 'to_date']) & pd.notnull(allo2.loc[:, 'from_date']))
    allo2.loc[idx, 'temp_days'] = (allo2.loc[idx, 'to_date'] - allo2.loc[idx, 'from_date'])
    allo2.loc[:, 'days'] = [e.days for e in allo2.loc[:, 'temp_days']]

    # the below appear to be an interal consents marker... and should not be included here as a replacement consent
    # is active at the same time at the consent with negitive number of days
    allo2.loc[allo2.loc[:, 'days'] < 0, 'days'] = 0

    allo2.loc[:, 'flux'] = allo2.loc[:, 'cav'] / 365 * allo2.loc[:, 'days'] / (end_time - start_time).days * -1
    allo2.loc[allo2.use_type == 'irrigation', 'flux'] *= 0.5

    out_data = allo2.reset_index().groupby('wap').aggregate({'flux': np.sum, 'crc': ','.join})
    out_data['consent'] = [tuple(e.split(',')) for e in out_data.loc[:, 'crc']]
    out_data = out_data.drop('crc', axis=1)
    out_data = out_data.dropna()

    out_data['type'] = 'well'
    out_data['zone'] = 's_wai'
    out_data.index.names = ['well']
    out_data.loc[:, 'flux'] *= 0.50  # start with a 50% scaling factor from CAV come back if time

    return out_data


def _check_waimak_wells():
    """
    some well checks
    :return:
    """
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp), index_col='crc')

    # option 2
    end_time = pd.datetime(2016, 12, 31)
    start_time = pd.datetime(2008, 1, 1)

    allo2 = allo.loc[np.in1d(allo['cwms'], ['Waimakariri']) &
                     (allo['take_type'] == 'Take Groundwater') &
                     ((allo.status_details.str.contains('Terminated')) | allo.status_details.str.contains('Issued'))]

    allo2.loc[:, 'to_date'] = pd.to_datetime(allo2.loc[:, 'to_date'], format='%d/%m/%Y', errors='coerce')
    allo2.loc[:, 'from_date'] = pd.to_datetime(allo2.loc[:, 'from_date'], format='%d/%m/%Y', errors='coerce')
    allo2.loc[allo2.loc[:, 'to_date'] > end_time, 'to_date'] = end_time
    allo2.loc[allo2.loc[:, 'to_date'] < start_time, 'to_date'] = None
    allo2.loc[allo2.loc[:, 'from_date'] < start_time, 'from_date'] = start_time
    allo2.loc[allo2.loc[:, 'from_date'] > end_time, 'from_date'] = None
    idx = (pd.notnull(allo2.loc[:, 'to_date']) & pd.notnull(allo2.loc[:, 'from_date']))
    allo2.loc[idx, 'temp_days'] = (allo2.loc[idx, 'to_date'] - allo2.loc[idx, 'from_date'])
    allo2.loc[:, 'days'] = [e.days for e in allo2.loc[:, 'temp_days']]

    # the below appear to be an interal consents marker... and should not be included here as a replacement consent
    # is active at the same time at the consent with negitive number of days
    allo2.loc[allo2.loc[:, 'days'] < 0, 'days'] = 0

    allo2.loc[:, 'flux'] = allo2.loc[:, 'cav'] / 365 * allo2.loc[:, 'days'] / (end_time - start_time).days * -1
    allo2.loc[allo2.use_type == 'irrigation', 'flux'] *= 0.5

    out_data = allo2.reset_index().groupby('wap').aggregate({'flux': np.sum, 'crc': ','.join})
    out_data['consent'] = [tuple(e.split(',')) for e in out_data.loc[:, 'crc']]
    out_data = out_data.drop('crc', axis=1)
    out_data = out_data.dropna()

    out_data['type'] = 'well'
    out_data['zone'] = 'n_wai'
    out_data.index.names = ['well']
    out_data.loc[:, 'flux'] *= 0.50  # start with a 50% scaling factor from CAV come back if time

    return out_data


def _get_s_wai_rivers():
    """
    get the well features that are being used to represent the selwyn hillfed streams
    :return:
    """
    # julian's report suggeststs that the hawkins is mostly a balanced flow and it was not included in
    # scott and thorley 2009 so it is not being simulated here

    # values for the upper selwyn, horoatat and waianiwanwa river are taken from scott and thorley 2009 and are evenly
    # distributed across the area as injection wells in layer 1

    # value in the shapefile reference the reach numbers from scott and thorley 2009

    no_flow = smt.get_no_flow(0)
    no_flow[no_flow < 0] = 0
    rivers = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/selwyn_hill_feds.shp".format(smt.sdp), 'reach', True)
    rivers[~no_flow.astype(bool)] = np.nan
    waian = pd.DataFrame(smt.model_where(rivers[np.newaxis, :, :] == 106), columns=['layer', 'row', 'col'])
    waian['well'] = ['waian{:04d}'.format(e) for e in waian.index]
    waian['flux'] = 0.142 * 86400 / len(waian)  # evenly distribute flux from scott and thorley 2009

    selwyn = pd.DataFrame(smt.model_where(rivers[np.newaxis, :, :] == 104), columns=['layer', 'row', 'col'])
    selwyn['well'] = ['selwyn{:04d}'.format(e) for e in selwyn.index]
    selwyn['flux'] = 4.152 * 86400 / len(selwyn)  # evenly distribute flux from scott and thorley 2009

    horo = pd.DataFrame(smt.model_where(rivers[np.newaxis, :, :] == 105), columns=['layer', 'row', 'col'])
    horo['well'] = ['horo{:04d}'.format(e) for e in horo.index]
    # evenly distribute flux from scott thorley 2009 but only including 7/32 of the river in my model
    horo['flux'] = 0.554 * 7. / 32 * 86400 / len(horo)

    # the hawkins 0.25 m3/s equally distributed between watsons bridge road and homebush road
    hawkins = pd.DataFrame(smt.model_where(rivers[np.newaxis, :, :] == 103), columns=['layer', 'row', 'col'])
    hawkins['well'] = ['hawkins{:04d}'.format(e) for e in hawkins.index]
    hawkins['flux'] = 0.25 * 86400 / len(hawkins)  # evenly distribute flux from scott and thorley 2009

    outdata = pd.concat((waian, selwyn, horo, hawkins))
    outdata['zone'] = 's_wai'
    outdata['type'] = 'river'
    outdata['consent'] = None

    return outdata


def get_s_wai_races():
    """
    get the wells that represent the race systems
    :return:
    """
    no_flow = smt.get_no_flow(0)
    race_array = smt.shape_file_to_model_array('{}/m_ex_bd_inputs/shp/s_wai_races.shp'.format(smt.sdp), 'race_code',
                                               True)
    race_array[np.isclose(no_flow, 0)] = np.nan
    nums = [1, 2, 3]
    names = ['ells_race', 'mal_race', 'papa_race']
    # for now the losses below are calculated as 80% loss over the races, and the influx volume is assumed to be the
    # baseflow from Selwyn District Council Water Race Management Plan 30 July 2013 table 2.1 the losses in the
    # ellesmere scheme is scaled by 38/156 because only a small portion of the scheme is in the model domain.
    min_flow = [1.539, 1.375, 1.231]  # excludes irrigation
    max_flow = [1.732, 2.210, 1.331]
    losses = [(i + j) / 2 for i, j in zip(min_flow, max_flow)]
    losses[0] *= 38 / 156
    proportion_lost = 0.89
    proportion_season = 0.75
    losses = [e * 86400 * proportion_lost * proportion_season for e in losses]
    outdata = pd.DataFrame()
    for num, name, loss in zip(nums, names, losses):
        idx = smt.model_where(np.isclose(race_array, num))
        keys = ['{}{:04d}'.format(name, e) for e in range(len(idx))]
        row = np.array(idx)[:, 0]
        col = np.array(idx)[:, 1]
        flux = loss / len(keys)
        temp = pd.DataFrame(index=keys, data={'row': row, 'col': col, 'flux': flux})
        outdata = pd.concat((outdata, temp))

    outdata['layer'] = 0
    outdata['zone'] = 's_wai'
    outdata['type'] = 'race'
    outdata['consent'] = None
    outdata.index.names = ['well']
    return outdata


def add_use_type(data):
    """
    add the use type to a set of well data
    :param data: well data (pd.DataFrame)
    :return:
    """
    data = deepcopy(data)
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    allo = allo.set_index('wap')
    data.loc[:, 'use_type'] = 'injection'
    for i in data.loc[data.type == 'well'].index:
        if np.in1d(np.atleast_1d(allo.loc[i, 'use_type']), ['stockwater', 'irrigation']).any():
            data.loc[i, 'use_type'] = 'irrigation-sw'
        else:
            data.loc[i, 'use_type'] = 'other'
    return data


if __name__ == '__main__':
    # tests
    from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.well_budget import get_well_budget
    from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.wells import get_max_rate, get_full_consent
    new = _get_wel_spd_v2()

    print('NEW')
    print(get_well_budget(new) / 86400)
    old = _get_wel_spd_v3()
    print('OLD')
    print(get_well_budget(old) / 86400)

    print('max_rate')
    max_rate = get_max_rate('opt')
    print(get_well_budget(max_rate)/86400)
    print('full CAV')
    full_cav = get_full_consent('opt')
    print(get_well_budget(full_cav)/86400)

    well_data = new
    well_data = well_data.loc[:, ['layer', 'row', 'col', 'flux', 'type']]
    well_data.to_csv(r"C:\Users\MattH\Desktop\to_brioch_2017_10_4/well_data.csv")

    new_nwai = _get_2014_2015_waimak_usage()
    raise
    nwai = get_nwai_wells()
    s_wells = _get_s_wai_wells()

    test = get_s_wai_races()
    old = _get_wel_spd_v1(recalc=False)
    new = _get_wel_spd_v3(recalc=True)

    n_wells_new = _check_waimak_wells()
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp), index_col='crc')
    chch_wells = _check_chch_wells()
    n_wells = get_nwai_wells()
    s_wells_all = _get_s_wai_wells()
    print 'done'
