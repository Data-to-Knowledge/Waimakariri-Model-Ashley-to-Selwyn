# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 21/07/2017 10:44 AM
"""
from __future__ import division
import env
import pandas as pd
import numpy as np
#from core.ecan_io import rd_sql, sql_db
from waimak_extended_boundry import smt
import geopandas as gpd
#import rasterio
import os
from future.builtins import input


#todo make this dataset avaible
def get_all_well_row_col(recalc=False):
    save_path = '{}/all_wells_row_col_layer2.csv'.format(smt.sdp)
    if os.path.exists(save_path) and not recalc:
        out_data = pd.read_csv(save_path,index_col=0)
        return out_data

    cont = input('are you sure you want to calculate all wells layer, row, col this takes several hours:\n {} \n continue y/n\n').lower()
    if cont != 'y':
        raise ValueError('script aborted')

    elv_sheet = pd.read_excel(
        env.sci('Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/xyz.xlsx'),
        index_Col=0)
    elv_sheet = elv_sheet.set_index('well')
    elv_sheet.loc[:,'accuracy_use'] = elv_sheet['ACCURACY (m)']

    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4) | (well_details_org['WMCRZone'] == 7) |
                                    (well_details_org['WMCRZone'] == 8)]  # keep only waimak selwyn and chch zones
    well_details = well_details[pd.notnull(well_details['DEPTH'])]
    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:, 'WELL_NO'] = [e.strip() for e in screen_details.loc[:, 'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    well_details.loc[:, 'WELL_NO'] = [e.strip() for e in well_details.loc[:, 'WELL_NO']]
    well_details = well_details.set_index('WELL_NO')

    rb = rasterio.open(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\dem.tif")

    out_data = pd.DataFrame(index=set(well_details.index))

    for well in out_data.index:
        out_data.loc[well, 'nztmx'] = well_details.loc[well, 'NZTMX']
        out_data.loc[well, 'nztmy'] = well_details.loc[well, 'NZTMY']
        out_data.loc[well, 'depth'] = well_details.loc[well, 'DEPTH']
        out_data.loc[well, 'num_screens'] = well_details.loc[well, 'Screens']
        out_data.loc[well,'aq_name'] = well_details.loc[well, 'AQUIFER_NAME']
        qarrl = out_data.loc[well, 'qar_rl'] = well_details.loc[well, 'QAR_RL']
        out_data.loc[well,'ref_ac'] = 0.1

        ref_level = well_details.loc[well, 'REFERENCE_RL']
        ground_ref_level = well_details.loc[well, 'GROUND_RL']
        if pd.isnull(ground_ref_level):
            ground_ref_level = 0
        if pd.isnull(ref_level):  # if there is no reference level assume it is at the ground from DEM
            if well in np.array(elv_sheet.index):
                ref_level = elv_sheet.loc[well, 'VALUE FOR USE']
                out_data.loc[well, 'ref_ac'] = elv_sheet.loc[well, 'accuracy_use']
            else:
                out_data.loc[well, 'ref_ac'] = 10
                x, y = well_details.loc[well, ['NZTMX', 'NZTMY']]
                ref_level = list(rb.sample([(x, y)]))[0][0] + ground_ref_level * -1
                if np.isclose(ref_level, -3.40282306e+38):
                    ref_level = np.nan
        out_data.loc[well, 'ref_level'] = ref_level
        out_data.loc[well, 'ground_level'] = ref_level + ground_ref_level

        if well in screen_details.index:
            top = np.min(screen_details.loc[well, 'TOP_SCREEN'])
            bot = np.max(screen_details.loc[well, 'BOTTOM_SCREEN'])
            out_data.loc[well, 'mid_screen_depth'] = (bot + top) / 2

            if out_data.loc[well,'num_screens'] == 2:
                top = np.atleast_1d(screen_details.loc[well, 'TOP_SCREEN']).max()
                bot = np.atleast_1d(screen_details.loc[well, 'BOTTOM_SCREEN']).min()
                out_data.loc[well,'distance_between_screen'] = top-bot

        else:
            out_data.loc[well, 'mid_screen_depth'] = well_details.loc[well, 'DEPTH'] - 2

    out_data.loc[:, 'mid_screen_elv'] = out_data.loc[:, 'ground_level'] - out_data.loc[:, 'mid_screen_depth']

    elv = smt.calc_elv_db()
    number_of_values = len(out_data.index)

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
    leapfrog_aq.loc[:,'use_aq_name'] = leapfrog_aq.loc[:,'aq_name']
    leapfrog_aq.loc[leapfrog_aq.use_aq_name.isnull(),'use_aq_name'] = leapfrog_aq.loc[leapfrog_aq.use_aq_name.isnull(),'aq_name_gi']

    for num,i in enumerate(out_data.index):
        if num%100 == 0:
            print ('completed {} of {}'.format(num,number_of_values))
        try:
            layer_by_depth, row, col = smt.convert_coords_to_matix(out_data.loc[i,'nztmx'],out_data.loc[i,'nztmy'],out_data.loc[i,'mid_screen_elv'])
            out_data.loc[i,'layer_by_depth'] = layer_by_depth
            out_data.loc[i,'row'] = row
            out_data.loc[i,'col'] = col

        except AssertionError as val:
            print(val)

    for well in out_data.index:
        try:
            out_data.loc[well, 'aquifer_in_confined'] = aq = leapfrog_aq.loc[well, 'use_aq_name']
            out_data.loc[well, 'layer_by_aq'] = aq_to_layer[aq]
        except KeyError:
            pass

    out_data.loc[:,'layer'] = out_data.loc[:, 'layer_by_depth']
    idx = out_data.layer_by_aq.notnull()
    out_data.loc[idx,'layer'] = out_data.loc[idx, 'layer_by_aq']

    for num, i in enumerate(out_data.index):
        row,col,layer = out_data.loc[i,['row', 'col', 'layer']]
        if any(pd.isnull([row,col,layer])):
            continue
        mx, my, mz = smt.convert_matrix_to_coords(int(row), int(col), int(layer), elv)
        out_data.loc[i, 'mx'] = mx
        out_data.loc[i, 'my'] = my
        out_data.loc[i, 'mz'] = mz

    out_data.to_csv(save_path)
    return out_data

if __name__ == '__main__':
    # note this takes some time to run
    all_well = get_all_well_row_col()
    print('done')