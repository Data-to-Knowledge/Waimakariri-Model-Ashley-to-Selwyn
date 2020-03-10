# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/07/2017 5:24 PM
"""

from __future__ import division
import env
from core.ecan_io import rd_sql, sql_db
import numpy as np
import pandas as pd
from core.classes.hydro import hydro
import rasterio
import geopandas as gpd
from copy import deepcopy
from waimak_extended_boundry import smt


def get_water_level_data(min_reading=1):
    elv_sheet = pd.read_excel(
        env.sci('Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/xyz.xlsx'),
        index_Col=0)
    elv_sheet = elv_sheet.set_index('well')

    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4) | (well_details_org['WMCRZone'] == 7) |
                                    (well_details_org['WMCRZone'] == 8)]  # keep only waimak selwyn and chch zones
    well_details = well_details[pd.notnull(well_details['DEPTH'])]

    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:, 'WELL_NO'] = [e.strip() for e in screen_details.loc[:, 'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    well_details.loc[:, 'WELL_NO'] = [e.strip() for e in well_details.loc[:, 'WELL_NO']]
    well_details = well_details.set_index('WELL_NO')

    data = hydro().get_data(mtypes=['gwl'], sites=list(well_details.index)).data
    data = data.loc['gwl'].reset_index()
    data1 = hydro().get_data(mtypes=['gwl_m'], sites=list(well_details.index)).data
    temp = data1.groupby(level=['mtype', 'site']).describe()[
        ['min', '25%', '50%', '75%', 'mean', 'max', 'count']].round(2)
    sites = list(temp.loc['gwl_m'].index[temp.loc['gwl_m']['count'] >= min_reading])

    data1 = data1.loc['gwl_m', sites].reset_index()
    data = pd.concat((data, data1), axis=0)
    data = data.set_index('site')
    data['month'] = [e.month for e in data.time]
    data['year'] = [e.year for e in data.time]

    rb = rasterio.open(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\dem.tif")

    out_data = pd.DataFrame(index=set(data.index))

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
                ref_level = elv_sheet.loc[well,'VALUE FOR USE']
                out_data.loc[well, 'ref_ac'] = elv_sheet.loc[well, 'ACCURACY (m)']
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
    out_data_org = deepcopy(out_data)
    data2008 = data.loc[data['year'] >= 2008]
    outputs = []
    for val, dat in zip(['_2008', '_all'], [data2008, data]):
        out_data = deepcopy(out_data_org)
        g = dat.loc[(dat.data > -999) & (dat.data < 999)].groupby('site')
        out_data['h2o_dpth_mean'] = g.aggregate({'data': np.mean})
        out_data['h2o_dpth_min'] = g.aggregate({'data': np.min})
        out_data['h2o_dpth_max'] = g.aggregate({'data': np.max})
        out_data['h2o_dpth_sd'] = g.aggregate({'data': np.std})
        out_data['readings_nondry'] = g.count().loc[:, 'data']

        for name, month in zip(['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'],
                               range(1, 13)):
            temp = dat.loc[(np.in1d(dat['month'], [month])) & (dat.data > -999) & (dat.data < 999)]
            tempg = temp.groupby('site')
            out_data['h2o_dpth_{}'.format(name)] = tempg.aggregate({'data': np.mean})
            out_data['reading_{}'.format(name)] = tempg.count().loc[:, 'data']

        g = dat.groupby('site')
        out_data['readings'] = g.count().loc[:, 'data']

        g = dat.loc[(dat.data <= -999)].groupby('site')
        out_data['count_dry'] = g.count().loc[:, 'data']

        g = dat.loc[(dat.data >= 999)].groupby('site')
        out_data['count_flowing'] = g.count().loc[:, 'data']

        keys_to_convert = pd.Series(out_data.keys())
        keys_to_convert = list(keys_to_convert[keys_to_convert.str.contains('h2o_dpth')])
        for key in keys_to_convert:
            out_data.loc[:, key.replace('dpth', 'elv')] = out_data.loc[:, key] + out_data.loc[:, 'ref_level']

        # add data quality flags
        all_readings = rd_sql(col_names=['WELL_NO', 'DATE_READ', 'TIDEDA_FLAG'],
                              where_col='WELL_NO',
                              where_val=[e.encode('ascii', 'ignore') for e in list(out_data.index)],
                              **sql_db.wells_db.daily)

        all_readings.loc[:, 'date'] = pd.to_datetime(all_readings.loc[:, 'DATE_READ'])
        if '2008' in val:
            all_readings = all_readings.loc[all_readings.date > pd.datetime(2008, 1, 1)]
        all_readings = all_readings.set_index('WELL_NO')

        # flag if more than 10% of the reading are by owners or more than 5% of the readings are suspect
        out_data.loc[:, 'owner_measured'] = False
        out_data.loc[:, 'bad_data'] = False

        for well in out_data.index:
            if well not in all_readings.index:
                continue
            flags = pd.Series(np.atleast_1d(all_readings.loc[well, 'TIDEDA_FLAG']))
            owner = flags.str.lower().str.contains('f').sum()
            bad_data = flags.str.lower().str.contains('x').sum()
            n_readings = len(flags)

            if owner / n_readings > 0.1:
                out_data.loc[well, 'owner_measured'] = True
            if bad_data / n_readings > 0.05:
                out_data.loc[:, 'bad_data'] = True

        # monitoring well need at least min of 6/year and not owner_measured or bad_data
        g = dat.groupby(['site', 'year'])
        temp = g.count()['data'].reset_index()
        temp = temp.loc[temp['data'] > 6]
        years = temp.groupby('site').count()['year']
        years = years.loc[years < 5]  # need at least 5 years of data to be considered a monitoring well

        out_data.loc[:, 'monitoring_well'] = False
        out_data.loc[list(years.index), 'monitoring_well'] = True
        out_data.loc[(out_data['bad_data']) | (out_data['owner_measured']), 'monitoring_well'] = False

        out_data.loc[:, 'water_range'] = out_data.loc[:, 'h2o_dpth_sep'] - out_data.loc[:, 'h2o_dpth_mar']

        out_data.loc[:, 'dry_percentage'] = (out_data.loc[:, 'readings'] - out_data.loc[:,
                                                                           'readings_nondry']) / out_data.loc[:,
                                                                                                 'readings'] * 100

        for name in ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']:
            out_data['s_factor_{}'.format(name)] = out_data['h2o_dpth_{}'.format(name)] / out_data.loc[:,
                                                                                          'h2o_dpth_mean']

        temp = out_data.loc[:, [u'reading_jan', u'reading_feb', u'reading_mar',
                                u'reading_apr', u'reading_may', u'reading_jun', u'reading_jul',
                                u'reading_aug', u'reading_sep', u'reading_oct', u'reading_nov',
                                u'reading_dec']]
        temp = temp.fillna(0)

        out_data.loc[:, 'samp_time_var'] = temp.std(axis=1) / temp.mean(axis=1)
        out_data = calc_target_offset(out_data)
        outputs.append(out_data)

    return outputs #2008 then all


def calc_target_offset(all_targets):
    vert_targets = pd.read_excel(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/Vertical gradient targets updated.xlsx"),
                                 sheetname='data_for_python', index_col=0)

    all_targets.loc[:,'vert_group'] = np.nan

    elv = smt.calc_elv_db()
    for i in vert_targets.index:
        all_targets.loc[i,'vert_group'] = vert_targets.loc[i,'group']
    number_of_values = len(all_targets.index)

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

    for num,i in enumerate(all_targets.index):
        if num%100 == 0:
            print ('completed {} of {}'.format(num,number_of_values))
        try:
            layer_by_depth, row, col = smt.convert_coords_to_matix(all_targets.loc[i,'nztmx'],all_targets.loc[i,'nztmy'],all_targets.loc[i,'mid_screen_elv'])
            all_targets.loc[i,'layer_by_depth'] = layer_by_depth
            all_targets.loc[i,'row'] = row
            all_targets.loc[i,'col'] = col

        except AssertionError as val:
            print(val)

    for well in all_targets.index:
        try:
            all_targets.loc[well, 'aquifer_in_confined'] = aq = leapfrog_aq.loc[well, 'use_aq_name']
            all_targets.loc[well, 'layer_by_aq'] = aq_to_layer[aq]
        except KeyError:
            pass

    all_targets.loc[:,'layer'] = all_targets.loc[:, 'layer_by_depth']
    idx = all_targets.layer_by_aq.notnull()
    all_targets.loc[idx,'layer'] = all_targets.loc[idx, 'layer_by_aq']

    for num, i in enumerate(all_targets.index):
        row,col,layer = all_targets.loc[i,['row', 'col', 'layer']]
        if any(pd.isnull([row,col,layer])):
            continue
        mx, my, mz = smt.convert_matrix_to_coords(int(row), int(col), int(layer), elv)
        all_targets.loc[i, 'mx'] = mx
        all_targets.loc[i, 'my'] = my
        all_targets.loc[i, 'mz'] = mz

    return all_targets

def define_error(outdata):
    # error as a percentage
    outdata.loc[:, 'total_error_m'] = np.nan
    outdata.loc[:, 'include_non-gap'] = False

    for well in outdata.index:
        readings = outdata.loc[well, 'readings_nondry']

        # add the different error terms:
        # measurement error
        me = np.abs(outdata.loc[well, 'h2o_dpth_sd'])
        if me==0:
            me = 0.01 * np.abs(outdata.loc[well, 'h2o_dpth_mean'])
        me = outdata.loc[well, 'm_error'] = 0.01 * np.abs(outdata.loc[well, 'h2o_dpth_mean'])
        # farmer vs ecan
        farm = outdata.loc[well, 'users_error'] = outdata.loc[well, 'owner_measured'] * 0.10 * np.abs(outdata.loc[
            well, 'h2o_dpth_mean'])
        # seasonal bias correction
        sea = outdata.loc[
            well, 'scorrection_error'] = outdata.loc[well, 'samp_time_var'] /3.3166 * 4.18
        # the above is based on the coefficient of variation(CV) of the monthly data/the cv if all sample are in 1 month
        #  * the mean plus one standard deviation of the water range (sept-mar) of all samples with as CV of less than 1
        # given time I would group these spatially, but I'm not that concerned. it will slightly over represent the
        # uncertainty of high seasonal CV in the coastal zone, and under-represent those with high seasonal CV in the
        # inland zone.  I plan to mostly uses only those with a seasonal CV of less than 0.1 for targets and gap fill
        # with others.

        # DEM error
        dem = outdata.loc[well, 'dem_error'] = outdata.loc[well,'ref_ac']

        low_rd_err = 0
        if readings < 20:
            low_rd_err = 15

        outdata.loc[well, 'low_rd_error'] = low_rd_err
        outdata.loc[well, 'total_error_m'] = dem + sea + me + (farm + low_rd_err) / (readings) ** 0.5

    idx = ((outdata.loc[:,'samp_time_var']<=0.1) &
           (outdata.loc[:, 'dem_error'] < 5) &
           (outdata.loc[:, 'readings'] > 20) &
           (outdata.loc[:,'total_error_m']<=2) &
           ((outdata.loc[:,'num_screens'] <= 1) | ((outdata.loc[:,'num_screens'] ==2) & (outdata.loc[:,'distance_between_screen']<5))))

    outdata.loc[idx,'include_non-gap'] = True

    return outdata

if __name__ == '__main__':
    path_2008 = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_woerror.csv")
    path_all = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_all_woerror.csv")
    recalc=True
    if recalc:
        outdata1, outdata2 = get_water_level_data()
        outdata1.to_csv(path_2008)
        outdata2.to_csv(path_all)
    else:
        outdata1 = pd.read_csv(path_2008)
        outdata2 = pd.read_csv(path_all)

    outpath_2008 = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_inc_error.csv")
    outpath_all = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_all_inc_error.csv")

    outdata1 = define_error(outdata1)
    outdata2 = define_error(outdata2)
    outdata1.to_csv(outpath_2008)
    outdata2.to_csv(outpath_all)


