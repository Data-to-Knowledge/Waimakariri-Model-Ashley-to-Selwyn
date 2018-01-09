# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/10/2017 2:04 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_at_wells import \
    get_hds_at_wells, _fill_df_with_bindata, _get_kstkpers, hds_no_data
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col
from core.ecan_io import sql_db, rd_sql
import os
import flopy
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.wells import \
    get_max_rate, get_full_consent
from copy import deepcopy
import rasterio


def _get_reliability_xyz(model_id):
    """
    returns a dataframe with index well list and x,y,z and i,j,k pump level, adiqute pen depth, use_pump_level
    z is the layer that use_pump_level comes from.
    :param model_id: which NSMC realisation
    :param well_list: the list of wells to use
    :return:
    """
    # get well list and location data
    well_list = pd.read_excel(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\Water supply wells\Well reliability assessment inputs.xlsx"))

    # add x,y / j,i
    all_wells = get_all_well_row_col()
    outdata = pd.merge(well_list, all_wells.loc[:, ['depth', 'i', 'j', 'nztmx', 'nztmy', 'ref_level', 'ground_level']],
                       how='right', right_index=True, left_index=True)

    # get pump level (top of screen 1, or if missing depth -2 m) # note need to convert to x,y,z
    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:, 'WELL_NO'] = [e.strip() for e in screen_details.loc[:, 'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    for well in outdata.index:
        if well in screen_details.index:
            top = np.min(screen_details.loc[well, 'TOP_SCREEN'])
        else:
            top = outdata.loc[well, 'depth'] - 2
        outdata.loc[well, 'pump_level'] = outdata.loc[well, 'ground_level'] - top

    # calculate addiquate penetration depth (ad_pen) #todo waiting on zeb
    outdata.loc[:, 'ad_pen'] = -99999  # todo a place holder, obviously

    # calculate use_pump_level
    outdata.loc[:, 'use_pump_level'] = outdata.loc[:, 'pump_level']
    idx = outdata.use_pump_level > outdata.ad_pen
    outdata[idx, 'use_pump_level'] = outdata.loc[idx, 'ad_pen']

    # add z/layer
    elv_db = smt.calc_elv_db()
    for well in outdata.index:
        row, col, elv = outdata.loc[well, ['i', 'j', 'use_pump_level']]
        outdata.loc[well, 'k'] = smt.convert_elv_to_k(row, col, elv, elv_db=elv_db)

    # get specific capacity data (specific_c)
    idx = (outdata.specific_c < 0.00210345) & (outdata.k ==0)
    layer0_sc = np.e**np.loadtxt(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_0.txt"))
    outdata.loc[idx,'specific_c'] = layer0_sc[outdata.loc[idx, 'i'], outdata.loc[idx, 'j']]

    idx = (outdata.specific_c < 0.00210345) & (outdata.k ==1)
    layer1_sc = np.e**np.loadtxt(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_1.txt"))
    outdata.loc[idx,'specific_c'] = layer1_sc[outdata.loc[idx, 'i'], outdata.loc[idx, 'j']]

    idx = (outdata.specific_c < 0.00210345) & (np.in1d(outdata.k, range(2,6)))
    layer2_5_sc = np.e**np.loadtxt(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_2-5.txt"))
    outdata.loc[idx,'specific_c'] = layer2_5_sc[outdata.loc[idx, 'i'], outdata.loc[idx, 'j']]

    idx = (outdata.specific_c < 0.00210345) & (np.in1d(outdata.k, range(6,11)))
    layer6_10_sc = np.e**np.loadtxt(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\inital_sc_data_rasters_extended\v2\arrays\ln_spe_capacity_layer_6-10.txt"))
    outdata.loc[idx,'specific_c'] = layer6_10_sc[outdata.loc[idx, 'i'], outdata.loc[idx, 'j']]

    # get pumping rate
    # 2.4 * cav rate (for 6 or 12 months) for irrigation takes and 1.2 * CAV rate (12 months) or max rate whichever is lower
    # for unconsented wells 10 m3/day
    increase_cav = get_full_consent(model_id)
    increase_cav.loc[
        increase_cav.use_type == 'irrigation-sw', 'flux'] *= 1.5 * 2  # peak month as 50% more water use and cav/6 months
    increase_cav.loc[increase_cav.use_type != 'irrigation-sw', 'flux'] *= 1.2
    max_rate = get_max_rate(model_id)
    idx = increase_cav.flux > max_rate.flux
    increase_cav.loc[idx, 'flux'] = max_rate.loc[idx, 'flux']
    outdata.loc[outdata.index, 'flux'] = increase_cav.loc[outdata.index, 'flux']

    return outdata


def get_model_well_reliability(model_id, model_path, indata):  # todo save this?
    """
    calculates the well reliability
    :param model_id: which NSMC realisation
    :param model_path: path to the model namefile with/without extension
    :param data: the outputs from _get_reliability_xyz which will speed up the for loops
    :return:
    """

    data = deepcopy(indata)
    # get simulation water level
    run_name = os.path.basename(model_path).replace('.nam', '').replace(model_id + '_',
                                                                        '')  # todo make sure this is used somewhere? or in next function up
    model_path = model_path.replace('.nam', '')
    hds_file = flopy.utils.HeadFile(model_path + '.hds')
    kstpkpers = _get_kstkpers(hds_file, kstpkpers=None, rel_kstpkpers=-1)
    kstpkper_names = 'model_water_level'
    data = _fill_df_with_bindata(hds_file, kstpkpers, kstpkper_names, data, hds_no_data, data)

    # adjust simulation to min water level #todo need zeb to provide more info

    # calculate drawdown and drawdown level #ignore the componenet of drawdown from average pumping in the cell, probably minor
    data.loc[:, 'dd_water_level'] = data.loc[:, 'low_water_level'] - ((data.loc[:, 'flux'] * 1000 / 86400) /
                                                                      data.loc[:, 'specific_c'])  # todo dbl check units
    # create reliability rating # todo handle if it is dry run by zeb
    temp = data.dd_water_level - data.use_pump_level
    data.loc[temp <= 0, 'rel_rate'] = 3
    data.loc[(temp > 0) & (temp <= 3), 'rel_rate'] = 2
    data.loc[temp > 3, 'rel_rate'] = 1

    data.loc[data.model_water_level < -777, 'rel_rate'] = 4  # already dry at average state #todo make sure hdry is set or compair to top/bottom this may not be needed if hdry is not set

    # create cost of reliability #todo handle if it is dry ask zeb how he wants to handle ask zeb what he wants as output
    # assign full consent volume.
    raise NotImplementedError
