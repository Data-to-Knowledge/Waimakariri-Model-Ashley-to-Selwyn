# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/10/2017 2:04 PM
"""

from __future__ import division
import env
import numpy as np
import pandas as pd
from glob import glob
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools.data_extraction.data_at_wells import \
    _fill_df_with_bindata, _get_kstkpers, hds_no_data
from waimak_extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col #todo move all wel col row to bc_data
from core.ecan_io import sql_db, rd_sql
import os
import flopy_mh as flopy
from waimak_extended_boundry.model_run_tools import \
    get_max_rate, get_full_consent
from copy import deepcopy


def _get_reliability_xyz(model_id, recalc=False):
    """
    returns a dataframe with index well list and x,y,z and i,j,k pump level, adiqute pen depth, use_pump_level
    z is the layer that use_pump_level comes from.
    :param model_id: which NSMC realisation
    :param recalc: the standard recalc to save in the temp pickle dir
    :return:
    """
    save_dir = os.path.join(smt.temp_pickle_dir, '{}_relability_xyz.hdf'.format(model_id))

    if os.path.exists(save_dir) and not recalc:
        outdata = pd.read_hdf(save_dir)
        return outdata

    # get well list and location data
    print('getting basic well data')
    well_list = pd.read_excel(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\Water supply wells\Well reliability assessment inputs.xlsx"),
        sheetname='Data for csv', index_col=0, names=['specific_c'])

    # add x,y / j,i
    all_wells = get_all_well_row_col()
    all_wells = all_wells.rename(columns={'row': 'i', 'col': 'j'})
    outdata = pd.merge(well_list, all_wells.loc[:, ['depth', 'i', 'j', 'nztmx', 'nztmy', 'ref_level', 'ground_level']],
                       how='left', right_index=True, left_index=True)

    # get pump level (top of screen 1, or if missing depth -2 m) # note need to convert to x,y,z
    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:, 'WELL_NO'] = [e.strip() for e in screen_details.loc[:, 'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    print('calculating pump level')
    for well in outdata.index:
        if well in screen_details.index:
            top = np.min(screen_details.loc[well, 'TOP_SCREEN'])
        else:
            top = outdata.loc[well, 'depth'] - 2
        outdata.loc[well, 'pump_level'] = outdata.loc[well, 'ground_level'] - top

    print('calculating adequate penetration')
    # calculate addiquate penetration depth (ad_pen) #
    # 2km wells drilled before 1 jan 2002, median depth then translate to elvation
    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.loc[(well_details.DATE_DRILLED.isnull() |
                                     (well_details.DATE_DRILLED < pd.datetime(2002, 1, 1))) &
                                    (well_details.WMCRZone == 4),  # keep only the waimakariri Zone
                                    ['NZTMX', 'NZTMY', 'DEPTH']]

    for well in outdata.index:
        temp = (
            (well_details.NZTMX - outdata.loc[well, 'nztmx']) ** 2 + (
                well_details.NZTMY - outdata.loc[well, 'nztmy']) ** 2)
        idx = temp <= 2000
        outdata.loc[well, 'ad_pen'] = outdata.loc[well, 'ground_level'] - np.nanmedian(well_details.loc[idx, 'DEPTH'])

    # calculate use_pump_level
    outdata.loc[:, 'use_pump_level'] = outdata.loc[:, 'pump_level']
    idx = outdata.use_pump_level > outdata.ad_pen
    outdata.loc[idx, 'use_pump_level'] = outdata.loc[idx, 'ad_pen']

    # add z/layer
    print('calculating z layer and getting specific capacity')
    elv_db = smt.calc_elv_db()
    for well in outdata.index:
        row, col, elv = outdata.loc[well, ['i', 'j', 'use_pump_level']]
        if pd.isnull(row) or pd.isnull(col) or pd.isnull(elv):
            continue
        outdata.loc[well, 'k'] = smt.convert_elv_to_k(int(row), int(col), elv, elv_db=elv_db)

        row, col, elv = outdata.loc[well, ['i', 'j', 'pump_level']]
        if pd.isnull(row) or pd.isnull(col) or pd.isnull(elv):
            continue
        outdata.loc[well, 'pump_k'] = smt.convert_elv_to_k(int(row), int(col), elv, elv_db=elv_db)

    # get specific capacity data (specific_c)
    outdata.loc[:,'sc_interpolated'] = False

    idx = (outdata.specific_c < 0.00210345) & (outdata.k == 0)
    layer0_sc = np.e ** np.loadtxt(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_0.txt"))
    outdata.loc[idx, 'specific_c'] = layer0_sc[outdata.loc[idx, 'i'].astype(int), outdata.loc[idx, 'j'].astype(int)]
    outdata.loc[idx, 'sc_interpolated'] = True

    idx = (outdata.specific_c < 0.00210345) & (outdata.k == 1)
    layer1_sc = np.e ** np.loadtxt(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_1.txt"))
    outdata.loc[idx, 'specific_c'] = layer1_sc[outdata.loc[idx, 'i'].astype(int), outdata.loc[idx, 'j'].astype(int)]
    outdata.loc[idx, 'sc_interpolated'] = True

    idx = (outdata.specific_c < 0.00210345) & (np.in1d(outdata.k, range(2, 6)))
    layer2_5_sc = np.e ** np.loadtxt(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_2-5.txt"))
    outdata.loc[idx, 'specific_c'] = layer2_5_sc[outdata.loc[idx, 'i'].astype(int), outdata.loc[idx, 'j'].astype(int)]
    outdata.loc[idx, 'sc_interpolated'] = True

    idx = (outdata.specific_c < 0.00210345) & (np.in1d(outdata.k, range(6, 11)))
    layer6_10_sc = np.e ** np.loadtxt(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_6-10.txt"))
    outdata.loc[idx, 'specific_c'] = layer6_10_sc[outdata.loc[idx, 'i'].astype(int), outdata.loc[idx, 'j'].astype(int)]
    outdata.loc[idx, 'sc_interpolated'] = True

    # get pumping rate
    print('getting pumping rate')
    # 2.4 * cav rate (for 6 or 12 months) for irrigation takes and 1.2 * CAV rate (12 months) or max rate whichever is lower
    # for unconsented wells 10 m3/day
    increase_cav = get_full_consent(model_id, missing_sd_wells=True)
    increase_cav.loc[
        increase_cav.use_type == 'irrigation-sw', 'flux'] *= 1  # peak month as 50% more water use and cav/6 months
    increase_cav.loc[increase_cav.use_type != 'irrigation-sw', 'flux'] *= 1
    max_rate = get_max_rate(model_id, missing_sd_wells=True)
    idx = increase_cav.flux < max_rate.flux
    increase_cav.loc[idx, 'flux'] = max_rate.loc[idx, 'flux']
    dups = increase_cav.index.duplicated()
    increase_cav = increase_cav.loc[dups]
    idx = set(outdata.index).intersection(increase_cav.index)
    outdata.loc[idx, 'flux'] = increase_cav.loc[idx, 'flux']

    outdata.loc[outdata.flux.isnull(), 'flux'] = -10  # assume wells not in increase_CAV are unconsented
    outdata.loc[:, 'flux'] *= -1  # make flux postive

    # add CAV
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(os.path.dirname(smt.sdp)))
    allo = allo.dropna(subset=['status_details'])
    allo = allo.loc[np.in1d(allo.status_details, ['Issued - Active', 'Issued - s124 Continuance']) & (
        allo.take_type == 'Take Groundwater')]
    missing_wells = pd.read_csv("{}/m_ex_bd_inputs/missing_sd_wells.csv".format(smt.sdp),
                                skiprows=3, index_col='well_no').loc[:, ['consent', 'max_rate', 'cav', 'use_type']]
    allo = allo.dropna(subset=['wap'])
    allo = allo.set_index('wap')
    allo = allo.loc[set(allo.index) - set(missing_wells.index)]
    allo = pd.concat((allo, missing_wells))
    allo = allo.reset_index()
    allo = allo.groupby('index').aggregate({'cav': np.sum})

    idx = set(outdata.index).intersection(allo.index)
    outdata.loc[idx, 'cav'] = allo.loc[idx, 'cav']
    outdata.loc[outdata.cav.isnull(), 'cav'] = 10 * 365  # assume that the CAV of missing wells is 10 m3/day * 365 days

    # for some reason there are wells which have a cav and not the flux
    idx = np.isclose(outdata.flux, 10) & ~np.isclose(outdata.cav, 3650)
    outdata.loc[idx, 'flux'] = 1.5 * outdata.loc[idx, 'cav'] / (365)  # assume these are poorly defined

    outdata = outdata.dropna(subset=[['k', 'i', 'j', 'pump_k']])
    # add cell bottom
    bots = elv_db[1:]
    outdata.loc[:, 'cell_bot'] = bots[[e for e in outdata.loc[:, ['k', 'i', 'j']].values.astype(int).transpose()]]
    outdata.loc[:, 'pump_cell_bot'] = bots[
        [e for e in outdata.loc[:, ['pump_k', 'i', 'j']].values.astype(int).transpose()]]

    outdata.to_hdf(save_dir, 'data', mode='w')

    return outdata


def get_model_well_reliability(model_path, indata):
    """
    calculates the well reliability
    :param model_path: path to the model namefile with/without extension
    :param indata: the outputs from _get_reliability_xyz which will speed up the for loops
    :return: run_name, data
    """

    data = deepcopy(indata)

    # get simulation water level
    run_name = os.path.basename(model_path).replace('.nam', '')
    model_path = model_path.replace('.nam', '')
    hds_file = flopy.utils.HeadFile(model_path + '.hds')
    kstpkpers = _get_kstkpers(hds_file, kstpkpers=None, rel_kstpkpers=-1)
    kstpkper_names = 'model_water_level'
    data = _fill_df_with_bindata(hds_file, kstpkpers, kstpkper_names, data, hds_no_data, data)
    assert not (data[kstpkper_names] < -777).any().any(), 'hdry must not be set for well reliablity calculations'

    # adjust simulation to min water level
    level_adj = smt.shape_file_to_model_array(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\Water supply wells\MeanvsMinWLadjustmentZones_v1.shp"),
        'WL_adj_m', False)
    data.loc[:, 'low_water_level'] = (data.loc[:, 'model_water_level'] +
                                      level_adj[[e for e in data.loc[:, ['i', 'j']].astype(int).values.transpose()]])
    # calculate drawdown and drawdown level
    # ignore the componenet of drawdown from average pumping in the cell, probably minor

    data.loc[:, 'dd_water_level'] = data.loc[:, 'low_water_level'] - ((data.loc[:, 'flux'] * 1000 / 86400) /
                                                                      data.loc[:, 'specific_c'])
    # create reliability rating
    temp = data.dd_water_level - data.use_pump_level
    data.loc[temp <= 0, 'ad_rel_rate'] = 3  # poor reliability
    data.loc[(temp > 0) & (temp <= 3), 'ad_rel_rate'] = 2  # moderate reliablity
    data.loc[temp > 3, 'ad_rel_rate'] = 1  # good reliability

    # completely unrelabile
    data.loc[data.model_water_level < data.cell_bot, 'ad_rel_rate'] = 4  # already dry at average state

    temp = data.dd_water_level - data.pump_level
    data.loc[temp <= 0, 'pump_rel_rate'] = 3  # poor reliability
    data.loc[(temp > 0) & (temp <= 3), 'pump_rel_rate'] = 2  # moderate reliablity
    data.loc[temp > 3, 'pump_rel_rate'] = 1  # good reliability

    # completely unrelabile
    data.loc[data.model_water_level < data.pump_cell_bot, 'pump_rel_rate'] = 4  # already dry at average state

    # create cost of reliability
    # assign full consent volume percentage of annual volume.  cost is in units of m3
    for inp, oup in zip(['ad_rel_rate', 'pump_rel_rate'], ['ad_cost', 'pump_cost']):
        for rel_rating, per in zip([1, 2, 3, 4], [0, .25, .5, 1]):
            idx = data[inp] == rel_rating
            data.loc[idx, oup] = data.cav * per

    data.loc[:, 'cutoff_el_rel_1'] = data.dd_water_level - 3.1
    data.loc[:, 'cutoff_el_rel_2'] = data.dd_water_level - 0.1
    data.loc[:, 'cutoff_el_rel_3'] = data.dd_water_level

    return run_name, data


def get_all_well_reliablity(indir, outdir):
    """
    get all well reliability for the models
    :param indir: the directory containing all of the models
    :param outdir: the directory to save the data
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with open(os.path.join(outdir, 'READ_ME.txt'), 'w') as f:
        f.write(
            """
            Below is an explanation of every variable in the datasheets in this folder
            
            specific_c: specific capactiy (l/s/m)
            depth: well depth (m)
            i: model row
            j: model column
            nztmx: longitude
            nztmy: latitude
            ref_level: well reference elevation (m)
            ground_level: well ground elevation (m)
            pump_level: pump elevation (top of screen) (m)
            ad_pen: adequate penetration elevation (m)
            use_pump_level: lower of ad_pen or pump_level elevation (m)
            k: model layer for use pump level
            pump_k: model layer pump level
            sc_interpolated: was SC interpolated
            flux: flux used m3/day
            cav: consented anual volume m3
            cell_bot: cell bottom for use_pump_level elevation (m)
            pump_cell_bot: bottom for pump_level elevation (m)
            model_water_level: elevation (m)
            low_water_level: model_water_level adjusted to expected annual lows elevation (m)
            dd_water_level: low_water_level + in well drawdown elevation (m)
            ad_rel_rate: reliability rating for adequate penetration scenario 1: good, 2: moderate, 3: poor, 4: dry at steady state
            pump_rel_rate:reliability rating for pump level only 1: good, 2: moderate, 3: poor, 4: dry at steady state
            ad_cost: cost for water volume in adiquate penetration scenario 
            pump_cost: cost for water volume in pumping senario
            cutoff_el_rel_1: cutoff for rating of 1 elevation (m)
            cutoff_el_rel_2: cutoff for rating of 2 elevation (m)
            cutoff_el_rel_3: cutoff for rating of 3 elevation (m)            
            """
        )

    model_paths = glob(os.path.join(indir, '*', '*.nam'))
    model_id = os.path.basename(model_paths[0]).split('_')[0]
    base_data = _get_reliability_xyz(model_id)
    # for loop run across
    for model_path in model_paths:
        name, outdata = get_model_well_reliability(model_path, base_data)
        outdata.to_csv(os.path.join(outdir, name + '.csv'))


if __name__ == '__main__':
    test = _get_reliability_xyz('NsmcBase',True)
    get_all_well_reliablity("D:\mh_waimak_models\NsmcBase_non_cc_forward_runs_2018-01-09",
                            "D:\mh_waimak_models\NsmcBase_non_cc_forward_runs_2018-01-09_well_reliablity")
