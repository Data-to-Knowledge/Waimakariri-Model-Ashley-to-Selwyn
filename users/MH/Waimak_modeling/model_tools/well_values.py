"""
Author: matth
Date Created: 26/04/2017 10:50 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import os
import pickle
from copy import deepcopy
from basic_tools import convert_coords_to_matix, calc_elv_db, convert_matrix_to_coords, get_well_postions, calc_z_loc
from users.MH.Waimak_modeling.supporting_data_path import sdp, temp_file_dir
from core.ecan_io import rd_sql,sql_db
from warnings import warn


def get_race_data(recalc = False):
    # load race data either from modflow or from pickle object
    pickle_path = "{}/inputs/pickeled_files/water_races.p".format(sdp)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
    else:
        elv_db = calc_elv_db()
        with open('{}/from_GNS/m_flooded/BH_OptMod_Flow.wel'.format(sdp)) as f:
            spd = f.readlines()[4:]
        spd = np.array([[float(j.strip()) for j in e.strip().split(' ')] for e in spd])
        spd[:,0:3] += -1
        outdata = pd.DataFrame(data=spd, columns=['layer', 'row', 'col', 'flux'])
        outdata = deepcopy(outdata[outdata.flux > 0]).reindex()
        outdata['type'] = 'race'
        names = ['race{:04d}'.format(e)for e in outdata.index]
        outdata['well'] = names
        outdata['consent'] = None

        for i in outdata.index:
            x,y,z = convert_matrix_to_coords(*outdata.loc[i,['row','col','layer']],
                                                                      elv_db=elv_db)
            outdata.loc[i, 'x'] = x
            outdata.loc[i, 'y'] = y
            outdata.loc[i, 'z'] = z

        # separate out boundry conditions
        bcs = pd.read_csv("{}/inputs/wells/boundry_fluxes.csv".format(sdp))  # this csv made from gis selection of races
        tempidx = np.in1d(outdata['well'],bcs['well'])

        outdata.loc[tempidx,'type'] = 'boundry_flux'
        pickle.dump(outdata,open(pickle_path,'w'))

    return outdata

def get_nwai_wells():
    # initalize and load extra data
    consent_details = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    consent_details2 = consent_details.set_index('crc')
    consent_details2 = consent_details2[(consent_details2['status_details']=='Issued - Active') | (consent_details2['status_details']=='Issued - s124 Continuance')]

    n_wai_wells = pd.read_csv("{}/inputs/wells/wells_Final_abstract_12_07.csv".format(sdp))  # data north of waimak (fouad)
    n_wai_wells = n_wai_wells.drop(['v1'],axis=1)
    n_wai_wells.columns = ['well','x', 'y', 'z', 'flux']
    n_wai_wells['type'] = 'well'
    n_wai_wells['consent'] = None

    for i in n_wai_wells.index:
        well = n_wai_wells.loc[i,'well']

        cons = tuple(consent_details2[consent_details2['wap']==well].index)

        n_wai_wells.set_value(i, 'consent', cons)
    return n_wai_wells

def _load_well_bc_data_from_csv_modflow():
    # initalize and load extra data
    consent_details = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    consent_details2 = consent_details.set_index('crc')
    consent_details2 = consent_details2[(consent_details2['status_details']=='Issued - Active') | (consent_details2['status_details']=='Issued - s124 Continuance')]
    outdata = pd.DataFrame(columns=['well','x', 'y', 'z', 'flux', 'type'])

    # load wells north of the waimakariri
    n_wai_wells = pd.read_csv("{}/inputs/wells/wells_Final_abstract_12_07.csv".format(sdp))  # data north of waimak (fouad)
    n_wai_wells = n_wai_wells.drop(['v1'],axis=1)
    n_wai_wells.columns = ['well','x', 'y', 'z', 'flux']
    n_wai_wells['type'] = 'well'
    n_wai_wells['consent'] = None

    for i in n_wai_wells.index:
        well = n_wai_wells.loc[i,'well']

        cons = tuple(consent_details2[consent_details2['wap']==well].index)

        n_wai_wells.set_value(i, 'consent', cons)

    outdata = pd.concat([outdata,n_wai_wells])

    # load wells south of the waimakariri
    s_wai_consents = pd.read_csv("{}/inputs/wells/consents_s_waimak.csv".format(sdp))  # consents south of waimak
    unknown = 0
    for i in s_wai_consents.index:
        con = s_wai_consents.loc[i,'consent']
        x, y, z, flux = s_wai_consents.loc[i,['x','y','z','flux']]
        wells = consent_details.wap[(consent_details.crc == con) & (consent_details.take_type=='Take Groundwater')]
        temp_data = pd.DataFrame(columns=['well','x', 'y', 'z', 'flux', 'type', 'consent'])
        temp_data['consent'] = temp_data['consent'].astype(object)


        if len(wells) == 0:
            temp_data.loc[0, ['well', 'x', 'y', 'z', 'flux', 'type']] = 'unknown{}'.format(unknown), x, y, z, flux, 'well'
            temp_data.set_value(0, 'consent', (con,))
            unknown += 1
            outdata = outdata.append(temp_data, ignore_index=True)
            continue


        for j, well in enumerate(wells):
            temp_data.loc[j,['well','x', 'y', 'z', 'flux', 'type']] = well, x, y, z, flux/len(wells), 'well'
            temp_data.set_value(j, 'consent', (con,))

        # remove any wells that are in nwai wells
        idx = np.in1d(temp_data.well,n_wai_wells.well,invert=True)
        outdata = outdata.append(temp_data[idx], ignore_index=True)

    # remove duplicates (ie well in 2 or more consents) and sum fluxes
    outdata['flux'] = outdata.groupby('well')['flux'].transform('sum')
    outdata = outdata.drop_duplicates('well')

    # add layer row col
    elv_db = calc_elv_db()
    for i in outdata.index:
        layer, row, col = convert_coords_to_matix(*outdata.loc[i,['x','y','z']], elv_db=elv_db)
        outdata.loc[i, 'layer'] = layer
        outdata.loc[i, 'row'] = row
        outdata.loc[i, 'col'] = col

    # add use and depth to these data
    usecodes = {'CI': 'com_industrial',
                'DA': 'dairy_use',
                'DO': 'domestic_supply',
                'DS': 'domestic_stockwater',
                'HP': 'heat_pump',
                'IB': 'sw_irrigation_bywash',
                'IR': 'irrigation',
                'OB': 'water_level_obs',
                'OT': 'other',
                'PG': 'power_gen',
                'PU': 'public_water_supply',
                'RE': 'recharge',
                'SC': 'small_community_supply',
                'ST': 'stock_supply',
                'TE': 'geo_investigation',
                'WQ': 'groundwater_quality',
                None: None,
                '': None}
    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details_org = well_details_org.set_index('WELL_NO')
    for i in outdata.index:
        well = outdata.loc[i,'well']
        if 'unknown' in well or 'race' in well:
            continue

        if well not in well_details_org.index:
            warn('{} not in well_details'.format(well))
            continue
        use1 = well_details_org.USE_CODE_1.str.strip().loc[well]
        use2 = well_details_org.USE_CODE_2.str.strip().loc[well]

        outdata.loc[i,'use1'] = usecodes[use1]
        outdata.loc[i,'use2'] = usecodes[use2]
        outdata.loc[i,'depth'] = well_details_org.loc[well,'DEPTH']


    # add water races
    race_data = get_race_data()
    outdata = outdata.append(race_data, ignore_index=True)

    outdata = outdata.set_index('well')
    return outdata

def get_original_well_data(recalc=False):
    """
    get the origiginal model wells as a dataframe with some additional parameters.  can either be calculated or read
    from a pickled file
    :param recalc: if True the prickled file will be re-calculated
    :return:
    """
    # load a pickled value or create a pickled value to load
    org_well_path = "{}/inputs/pickeled_files/original_well_data.p".format(sdp)
    if (os.path.exists(org_well_path)) and (not recalc):
        well_data = pickle.load(open(org_well_path))
    else:
        well_data = _load_well_bc_data_from_csv_modflow()
        pickle.dump(well_data,open(org_well_path,'w'))

    return well_data

def convert_well_data_to_stresspd(well_data_in):
    # convert a dataframe of well features (x,y,z,flux,type, etc.) to well standard stress period data
    # do something similar for concentration data?

    # do groupby statistics to sum everything that is in the same layer, col, row
    g = well_data_in.groupby(['layer','row','col'])
    well_ag = g.aggregate({'flux': np.sum}).reset_index()

    outdata = []
    for i in well_ag.index:
        outdata.append(list(well_ag.loc[i,['layer','row','col','flux']]))

    return outdata


def get_all_well_data(recalc=False):
    # return a dataframe of useful information for all of the wells in the waimak, selwyn, chch
    pickle_path = "{}/inputs/pickeled_files/all_well_indexs.p".format(sdp)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata


    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4) | (well_details_org['WMCRZone'] == 7) |
                                (well_details_org['WMCRZone'] == 8)]  # keep only waimak selwyn and chch zones
    well_details = well_details[well_details['Well_Status'] == 'AE']
    well_details = well_details[pd.notnull(well_details['DEPTH'])]

    # load all wells in the model
    outdata = get_original_well_data()

    # create new data frame
    outdata2 = pd.DataFrame(columns=['well', 'x', 'y', 'z', 'layer', 'row', 'col', 'flux', 'type'])

    wells_to_add = well_details['WELL_NO'][np.in1d(well_details['WELL_NO'],outdata.index,invert=True)]
    outdata2['well'] = wells_to_add
    outdata2['flux'] = np.nan
    outdata2['type'] = 'not_added'

    elp = '{}/well_position_errors'.format(temp_file_dir)
    # calculate layer,row,col
    layer, row, col = get_well_postions(wells_to_add,screen_handling='middle', raise_exct=False,
                                        error_log_path=elp, one_val_per_well=True)
    outdata2['layer'] = layer
    outdata2['row'] = row
    outdata2['col'] = col

    outdata2 = outdata2[pd.notnull(outdata2['row'])]

    elv_db = calc_elv_db()
    for i in outdata2.index:
        x, y, z = convert_matrix_to_coords(*outdata2.loc[i,['row','col','layer']],elv_db=elv_db)
        outdata2.loc[i,'x'] = x
        outdata2.loc[i,'y'] = y
        outdata2.loc[i,'z'] = z

    outdata2 = outdata2.set_index('well')

    outdata_final = pd.concat([outdata, outdata2])

    #add user codes
    outdata_final['use1'] = None
    outdata_final['use2'] = None
    outdata_final['well_type'] = None
    outdata_final['screen_z'] = None
    outdata_final['depth'] = None
    usecodes = {'CI': 'com_industrial',
                'DA': 'dairy_use',
                'DO': 'domestic_supply',
                'DS': 'domestic_stockwater',
                'HP': 'heat_pump',
                'IB': 'sw_irrigation_bywash',
                'IR': 'irrigation',
                'OB': 'water_level_obs',
                'OT': 'other',
                'PG': 'power_gen',
                'PU': 'public_water_supply',
                'RE': 'recharge',
                'SC': 'small_community_supply',
                'ST': 'stock_supply',
                'TE': 'geo_investigation',
                'WQ': 'groundwater_quality',
                None: None,
                '': None}

    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    well_details_org = well_details_org.set_index("WELL_NO")
    for well in outdata_final.index:
        if 'unknown' in well or 'race' in well:
            continue

        if well not in well_details_org.index:
            warn('{} not in well_details'.format(well))
            continue
        use1 = well_details_org.USE_CODE_1.str.strip().loc[well]
        use2 = well_details_org.USE_CODE_2.str.strip().loc[well]

        outdata_final.loc[well,'use1'] = usecodes[use1]
        outdata_final.loc[well,'use2'] = usecodes[use2]
        outdata_final.loc[well,'depth'] = well_details_org.loc[well,'DEPTH']
        outdata_final.loc[well, 'well_type'] = well_details_org.WELL_TYPE.loc[well]

        # calculate z screens
        row = int(outdata_final.loc[well, 'row'])
        col = int(outdata_final.loc[well, 'col'])

        # get screen elevations
        screen_num = well_details_org.loc[well, 'Screens']
        ref_level = well_details_org.loc[well, 'REFERENCE_RL']
        ground_level = well_details_org.loc[well, 'GROUND_RL'] + ref_level
        if pd.isnull(ground_level):
            ground_level = elv_db[0, row, col]  # take as top of model

        if screen_num == 0:
            screen_elvs = [outdata_final.loc[well,'z']]
        else:
            screen_elvs = []
            for j in range(1, screen_num + 1):
                top = ground_level - screen_details['TOP_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                   (screen_details['SCREEN_NO'] == j)].iloc[0]
                bottom = ground_level - screen_details['BOTTOM_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                         (screen_details['SCREEN_NO'] == j)].iloc[0]
                middle = ground_level - (top + bottom) / 2
                screen_elvs.append(top)
                screen_elvs.append(middle)
                screen_elvs.append(bottom)

        outdata_final = outdata_final.set_value(well,'screen_z',tuple(screen_elvs))



    pickle.dump(outdata_final,open(pickle_path,'w'))

    return outdata_final

def get_model_well_full_consented(recalc=False):
    """
    replaces the flux value with the full consented volumes for the wells north of the waimakariri
    :param recalc:
    :return:
    """
    pickle_path = "{}/inputs/pickeled_files/model_well_consented_vol.p".format(sdp)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    outdata = get_original_well_data()
    consent_data = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    consent_data = consent_data.set_index('crc')
    consent_data = consent_data.drop_duplicates()
    for well in get_n_wai_well_list():
        con = outdata.loc[well,'consent']
        if len(con) == 0:
            outdata.loc[well, 'flux'] = 0
            continue

        cdata = consent_data.loc[con,['cav','wap']]
        cdata = cdata['cav'][cdata['wap']==well].sum()
        outdata.loc[well,'flux'] = cdata/-365

    pickle.dump(outdata, open(pickle_path, 'w'))

    return outdata

def get_model_well_max_rate(recalc=False):
    """
    replaces the flux value with the full consented volumes for the wells north of the waimakariri
    :param recalc:
    :return:
    """
    pickle_path = "{}/inputs/pickeled_files/model_well_max_rate.p".format(sdp)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    outdata = get_original_well_data()
    consent_data = pd.read_csv("{}/inputs/wells/allo_gis.csv".format(sdp))
    consent_data = consent_data.set_index('crc')
    consent_data = consent_data.drop_duplicates()
    for well in get_n_wai_well_list():
        con = outdata.loc[well,'consent']
        if len(con) ==0:
            outdata.loc[well, 'flux'] = 0
            continue

        # first try to pull the max rate for the wap
        cdata = consent_data.loc[con,['max_rate_wap','wap']]
        cdata = cdata['max_rate_wap'][cdata['wap']==well].sum()

        if pd.isnull(cdata):  # if the max rate for the wap is not defined take the pro-rataed max_rate
            cdata = consent_data.loc[con,['max_rate','wap']]
            cdata = cdata['max_rate'][cdata['wap']==well].sum()

        if pd.isnull(cdata):
            cdata = consent_data.loc[con,['cav','wap']]
            cdata = cdata['cav'][cdata['wap']==well].sum()
            outdata.loc[well,'flux'] = cdata/-150
        else:
            outdata.loc[well, 'flux'] = cdata*-86.4

    pickle.dump(outdata, open(pickle_path, 'w'))

    return outdata

def get_n_wai_well_list():
    n_wai_wells = pd.read_csv("{}/inputs/wells/wells_Final_abstract_12_07.csv".format(sdp))  # data north of waimak (fouad)
    n_wai_wells = n_wai_wells.drop(['v1'],axis=1)
    n_wai_wells.columns = ['well','x', 'y', 'z', 'flux']

    return list(set(n_wai_wells.loc[:,'well']))


if __name__ == '__main__':
    from osgeo import gdal
    test = get_race_data(True)
    #temp = get_original_well_data(True)
    #temp2 = convert_well_data_to_stresspd(temp)

    #temp3 = get_all_well_data(True)
    #temp4 = get_model_well_full_consented(True)
    #temp5 = get_model_well_max_rate(True)
    #print 'done'

