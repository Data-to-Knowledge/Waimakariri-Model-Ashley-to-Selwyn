# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/09/2017 3:55 PM
"""

from __future__ import division

import os
import pickle

import numpy as np
import pandas as pd

from core.ecan_io import rd_sql, sql_db
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import _get_wel_spd_v1, _get_wel_spd_v2
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.LSR_arrays import \
    get_ird_base_array, get_lsr_base_period_inputs
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_base_well, temp_pickle_dir
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.well_budget import get_well_budget
from users.MH.Waimak_modeling.supporting_data_path import sdp
from copy import deepcopy
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col

_missing_well_use_converter = {'Aquaculture': 'other',
                               'Arable Farming': 'irrigation-sw',
                               'Dewatering': 'other',
                               'Intensive Farming': 'irrigation-sw',
                               'Irrigation': 'irrigation-sw',
                               'Irrigation and community supply': 'irrigation-sw',
                               'Irrigation, domestic and Stockwater': 'irrigation-sw',
                               'Other': 'other',
                               'Public Water Supply (Municipal/Community)': 'other',
                               'Recreation': 'other',
                               'Stockwater': 'irrigation-sw'}

# for stream depletion things
def get_race_data(model_id):
    """
    all influx wells in the well data (e.g. streams, races, boundary fluxes)
    :param model_id: the NSMC realisation
    :return:
    """
    outdata = get_base_well(model_id, True)  # use original pumping wells becasue it doesn't matter
    outdata = outdata.loc[outdata.loc[:, 'type'] != 'well']
    return outdata


def get_full_consent(model_id, org_pumping_wells=False, missing_sd_wells=False,
                     recalc=False):
    """
    EAV/consented annual volume I need to think about this one
    :param model_id: the NSMC realisation
    :param org_pumping_wells: Boolean if true use the model period wells else use the 2014/2015 pumping
    :param missing_sd_wells: bool if True include the wells that mark found missing
    :param recalc: Boolean, if true recalculate otherwise load from a pickle
    :return: the full dataframe with just the flux converted
    """
    if org_pumping_wells:
        new = 'mod_period'
    else:
        new = '2014_15'
    missing_sd = '_inc_missing_sd' if missing_sd_wells else ''

    pickle_path = "{}/model_{}_well_full_abstraction_{}{}.p".format(temp_pickle_dir, model_id, new, missing_sd)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    outdata = get_base_well(model_id, org_pumping_wells=org_pumping_wells)
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    allo = allo.dropna(subset=['status_details'])
    allo = allo.loc[np.in1d(allo.status_details, ['Issued - Active', 'Issued - s124 Continuance']) & (
        allo.take_type == 'Take Groundwater')]

    # adjust irrigation wells cav where they are calcualted on the July to June year to an October to April?
    idx = (allo.use_type == 'irrigation') & (allo.from_month != 'OCT')  # nans assumed to be july to june
    allo.loc[idx, 'cav'] *= 6 / 12
    fluxes = allo.groupby('wap').aggregate({'cav': np.sum})
    allo2 = allo.drop_duplicates('wap')
    allo2 = allo2.set_index('wap')
    fluxes = pd.merge(fluxes, pd.DataFrame(allo2.loc[:, 'cwms']), right_index=True, left_index=True)

    print_digong = False  # something I did to debug, it looks like there are some differences in teh well list,
    # but they look ok
    if print_digong:
        print ('only waimak')
        print(len(set(fluxes.loc[np.in1d(fluxes.cwms, ['Waimakariri'])].index) - set(outdata.index)))
        print((set(fluxes.loc[np.in1d(fluxes.cwms, ['Waimakariri'])].index) - set(outdata.index)))
        print ('only chch')
        print(len(set(fluxes.loc[np.in1d(fluxes.cwms, ['Christchurch - West Melton'])].index) - set(outdata.index)))
        print((set(fluxes.loc[np.in1d(fluxes.cwms, ['Christchurch - West Melton'])].index) - set(outdata.index)))
        print('number of missing in model')
        print(len(set(outdata.loc[outdata.type == 'well'].index) - set(fluxes.index)))

    outdata.loc[outdata.type == 'well', 'flux'] = 0
    idx = set(fluxes.index).intersection(outdata.index)
    outdata.loc[idx, 'flux'] = fluxes.loc[idx, 'cav'] / -365

    if missing_sd_wells:
        # there are a couple of duplicates e.g 1 well on two consents
        # these will be handled when I amalgamate up to the model grid
        all_wells = get_all_well_row_col()
        idx = pd.notnull(all_wells.loc[:, ['layer', 'row', 'col']]).all(axis=1)
        missing_wells = pd.read_csv("{}/m_ex_bd_inputs/missing_sd_wells.csv".format(smt.sdp),
                                    skiprows=3, index_col='well_no').loc[:, ['consent', 'max_rate', 'cav', 'use_type']]
        missing_wells.loc[missing_wells.use_type.isnull(), ['use_type']] = 'other'
        missing_wells = missing_wells.replace({'use_type':_missing_well_use_converter})

        missing_wells['max_rate'] *= -86.4  # to convert from l/s to m3/day
        missing_wells = pd.merge(missing_wells, all_wells.loc[idx, ['layer', 'row', 'col']],
                                 right_index=True, left_index=True)
        missing_wells.loc[:, 'type'] = 'well'
        missing_wells.loc[:, 'zone'] = 'n_wai'
        missing_wells.loc[idx, 'flux'] = missing_wells.loc[idx, 'cav'] / -365
        idx = ~np.in1d(outdata.index, missing_wells.index)
        outdata = pd.concat((outdata.loc[idx], missing_wells))

    pickle.dump(outdata, open(pickle_path, 'w'))
    return outdata


def get_max_rate(model_id, org_pumping_wells=False, missing_sd_wells=False,
                 recalc=False):
    """
    replaces the flux value with the full consented volumes for the wells north of the waimakariri
    :param model_id: the NSMC realisation
    :param org_pumping_wells: Boolean if true use the model period wells else use the 2014/2015 pumping
    :param missing_sd_wells: bool if True include the wells that mark flagged as missing from the waimakariri zone
    :param recalc: Boolean, if true recalculate otherwise load from a pickle
    :return: the full dataframe with just the flux converted
    """
    if org_pumping_wells:
        new = 'mod_period'
    else:
        new = '2014_15'

    missing_sd = '_inc_missing_sd' if missing_sd_wells else ''

    pickle_path = "{}/model_{}_well_max_rate_{}{}.p".format(temp_pickle_dir, model_id, new, missing_sd)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    outdata = get_base_well(model_id, org_pumping_wells)
    allo = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    allo = allo.dropna(subset=['status_details'])
    allo = allo.loc[np.in1d(allo.status_details, ['Issued - Active', 'Issued - s124 Continuance']) & (
        allo.take_type == 'Take Groundwater')]

    # adjust irrigation wells cav where they are calcualted on the July to June year to an October to April?
    idx = (allo.use_type == 'irrigation') & (allo.from_month != 'OCT')  # nans assumed to be july to june
    allo.loc[idx, 'cav'] *= 6 / 12
    fluxes = allo.groupby('wap').aggregate({'cav': np.sum, 'max_rate_wap': np.sum, 'max_rate': np.sum})
    allo2 = allo.drop_duplicates('wap')
    allo2 = allo2.set_index('wap')
    fluxes = pd.merge(fluxes, pd.DataFrame(allo2.loc[:, 'cwms']), right_index=True, left_index=True)

    # first max_rate_wap then max_rate then cav/365
    fluxes.loc[:, 'flux'] = fluxes.loc[:, 'max_rate_wap'] * -86.4  # to convert from l/s to m3/day
    fluxes.loc[fluxes.flux.isnull(), 'flux'] = fluxes.loc[:, 'max_rate'] * -86.4
    fluxes.loc[fluxes.flux.isnull(), 'flux'] = fluxes.loc[:, 'cav'] / -365  # this could cause an inverse of max rate
                                                                            # and 150 day rate
    outdata.loc[outdata.type == 'well', 'flux'] = 0
    idx = set(fluxes.index).intersection(outdata.index)
    outdata.loc[idx, 'flux'] = fluxes.loc[idx, 'flux']

    if missing_sd_wells:
        # there are a couple of duplicates e.g 1 well on two consents
        # these will be handled when I amalgamate up to the model grid

        all_wells = get_all_well_row_col()
        idx = pd.notnull(all_wells.loc[:, ['layer', 'row', 'col']]).all(axis=1)
        missing_wells = pd.read_csv("{}/m_ex_bd_inputs/missing_sd_wells.csv".format(smt.sdp),
                                    skiprows=3, index_col='well_no').loc[:, ['consent', 'max_rate', 'cav', 'use_type']]
        missing_wells.loc[missing_wells.use_type.isnull(), ['use_type']] = 'other'
        missing_wells = missing_wells.replace({'use_type':_missing_well_use_converter})

        missing_wells['flux'] = missing_wells.loc[:, 'max_rate'] * -86.4  # to convert from l/s to m3/day
        missing_wells = pd.merge(missing_wells, all_wells.loc[idx, ['layer', 'row', 'col']],
                                 right_index=True, left_index=True)
        missing_wells.loc[:, 'type'] = 'well'
        missing_wells.loc[:, 'zone'] = 'n_wai'
        idx = ~np.in1d(outdata.index, missing_wells.index)
        outdata = pd.concat((outdata.loc[idx], missing_wells))

    pickle.dump(outdata, open(pickle_path, 'w'))
    return outdata


def get_forward_wells(model_id, full_abstraction=False, cc_inputs=None, naturalised=False, full_allo=False,
                      pc5=False, org_pumping_wells=False):
    """
    gets the pumping data for the forward runs
    :param model_id: which NSMC realisation to use
    :param full_abstraction: boolean use the CAV (think about what happens with irrigation abstraction)
    :param cc_inputs: use these to apply scaling factors for the pumping (think about how to work with these spatially)
    :param naturalised: boolean, if True use only the fixed inputs (e.g. rivers, boundary fluxes.  No races)
    :param full_allo: boolean, if True scale the wells by the amount allocated in each zone (could be a dictionary of boolean for each subzone)
    :param org_pumping_wells: if True use the model peiod wells if false use the 2014-2015 usage
    :return:
    """
    new_water_needed = 0
    # check input make sense
    if full_abstraction and naturalised:
        raise ValueError('cannot both fully abstracted and naturalised')
    if full_allo and naturalised:
        raise ValueError('cannot both be fully allocated and naturalised')

    outdata = get_base_well(model_id, org_pumping_wells)
    if not org_pumping_wells:
        max_pumping = get_full_consent(model_id, org_pumping_wells)
        idx = outdata.index
        outdata.loc[(outdata.loc[idx, 'flux'] < max_pumping.loc[idx, 'flux']) & (outdata.cwms == 'waimak') & (
            outdata.type == 'well'), 'flux'] = max_pumping.loc[idx, 'flux']

    if full_abstraction:
        idx = outdata.loc[(outdata.type == 'well') & (outdata.cwms == 'waimak')].index
        outdata.loc[idx, 'flux'] = get_full_consent(model_id, org_pumping_wells).loc[idx, 'flux']
    else:
        if pc5 and not full_abstraction:
            outdata.loc[(outdata.loc[:, 'use_type'] == 'irrigation-sw') & (outdata.cwms == 'waimak'), 'flux'] *= 3 / 4
            # an inital 1/4 reduction for pc5 to
            # account for the decreased irrgation demand for with more efficent irrigation this number comes from
            # prorataing the difference between 80% and 100% irrigation LSRM outputs to the percentage of irrigation
            # from sw and gw in the zones

    if full_allo:
        allo_mult = get_full_allo_multipler(org_pumping_wells)
        idx = allo_mult.index
        outdata.loc[idx, 'flux'] *= allo_mult.loc[idx]

    if naturalised:
        outdata = outdata.loc[
            np.in1d(outdata.loc[:, 'type'], ['boundry_flux', 'llr_boundry_flux', 'river', 'ulr_boundry_flux'])]
    cc_mult = 1
    if cc_inputs is not None:
        temp = deepcopy(cc_inputs)
        try:
            temp.pop('cc_to_waimak_only')
        except KeyError:
            pass
        if all(pd.isnull(temp.values())):
            pass
        elif any(pd.isnull(temp.values())):
            raise ValueError('null and non-null values returned for cc_inputs')
        else:
            if org_pumping_wells:
                raise NotImplementedError('model period wells not implemented for cc senarios due to climate change '
                                          'problems around the multiplier and the CAV')

            cc_mult = get_cc_pumping_muliplier(cc_inputs)  # only apply cc multiplier to the waimakariri zone
            outdata.loc[(outdata.loc[:, 'use_type'] == 'irrigation-sw') & (outdata.cwms == 'waimak'), 'flux'] *= cc_mult
            temp = outdata.loc[:, 'flux'].sum()
            # pumping is truncated at full allocation and abstraction value but only for the waimakariri zone, selwyn/chch noth changed
            # we assume that any additional irrigation demand would be met with surface water schemes from the alpine rivers

            max_pumping = get_full_consent(model_id, org_pumping_wells)
            allo_mult = get_full_allo_multipler(org_pumping_wells)
            idx = allo_mult.index
            max_pumping.loc[idx, 'flux'] *= allo_mult.loc[idx]

            idx = outdata.index
            # less than as the fluxes are negative
            outdata.loc[(outdata.loc[idx, 'flux'] < max_pumping.loc[idx, 'flux']) & (outdata.cwms == 'waimak') &
                        (outdata.type == 'well'), 'flux'] = max_pumping.loc[idx, 'flux']
            # but may not apply as we'll only do the forward CC runs with 2014/15 pumping
            new_water_needed = (temp - outdata.loc[:, 'flux'].sum()) * -1

    return outdata, cc_mult, new_water_needed


def get_cc_pumping_muliplier(cc_inputs):
    """
    get the relative pumping multipler from the irrigation demand only calculated over the waimakariri zone
    :param cc_inputs: the cc inputs dictionary see forward wells
    :return: float
    """
    # return a single value for now which is senario/baseline for the senario see below
    amalg_dict = {None: 'mean', 'mean': 'mean', 'tym': 'period_mean', 'low_3_m': '3_lowest_con_mean',
                  'min': 'lowest_year'}

    sen = 'current'  # irrigation demand does not change with irrigation efficiency
    rcp = cc_inputs['rcp']
    rcm = cc_inputs['rcm']
    per = cc_inputs['period']
    at = amalg_dict[cc_inputs['amag_type']]
    # base period for pumping multiplier is RCPpast period mean for the same senario and rcm
    ird_current_period = get_ird_base_array(*get_lsr_base_period_inputs(sen, rcp, rcm, per, at))
    ird_modeled_period = get_ird_base_array(sen, rcp, rcm, per, at)
    outdata = ird_modeled_period / ird_current_period
    w_idx = get_zone_array_index(
        'waimak')  # regardless where the cc input are applied only apply ccmul to waimaik so only use that zone
    outdata = outdata[w_idx]
    return np.nanmean(outdata)


def get_full_allo_multipler(org_pumping_wells, recalc=False):
    """
    get the multipliers (on a per well basis) to move to full allocation
    :param org_pumping_wells: if true use model period wells, else use 2014/15 wells
    :param recalc: usual recalc
    :return: pd.Series
    """
    # return a series with index well numbers and values multiplier
    # these wells only in Waimakariri Zone
    # can get this from gw allocation zone in well_details
    """ Eyre	99.1	100.5	101% # JUST CONSIDER THIS 1
        Cust	56.3	15.5	27%
        Ashley	29.4	15.5	53%
        FROM WAIMAKARIRI CURRENT STATE REPORT
        https://punakorero/groups/plansec/WaimakAsh/research/Current%20State
        /Current%20State%20Groundwater%20Quantity%20Waimakariri_DRAFT_RevC.docx?web=1
    """
    if org_pumping_wells:
        new = 'mod_period'
    else:
        new = '2014_15'

    pickle_path = "{}/model_well_allo_mult_{}.p".format(smt.pickle_dir, new)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.set_index('WELL_NO')
    if org_pumping_wells:
        wells = _get_wel_spd_v1()
    else:
        wells = _get_wel_spd_v2()
    wells = pd.merge(wells, pd.DataFrame(well_details.loc[:, 'AllocationZone']), right_index=True, left_index=True)
    wells = wells.rename(columns={'AllocationZone': 'a_zone'})
    allo_zones = rd_sql(**sql_db.wells_db.allo_zones)
    replacement_dict = dict(zip(allo_zones.ZoneCode, allo_zones.ZONE))
    wells = wells.replace({'a_zone': replacement_dict})
    wells.loc[:, 'mult'] = 1
    wells.loc[wells.a_zone == 'Ashley', 'mult'] = 100 / 53  # see above for citation for these numbers
    wells.loc[wells.a_zone == 'Cust', 'mult'] = 100 / 27
    outdata = wells.loc[:, 'mult']

    pickle.dump(outdata, open(pickle_path, 'w'))
    return outdata


if __name__ == '__main__':
    wells = get_full_consent('NsmcBase', missing_sd_wells=True,recalc=True)
    wells = get_max_rate('NsmcBase', missing_sd_wells=True,recalc=True)
    print('done')
