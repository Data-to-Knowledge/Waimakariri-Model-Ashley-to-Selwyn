"""
coding=utf-8
Author: matth
Date Created: 26/06/2017 2:35 PM

depreciated due to chahagnes in core
"""

from __future__ import division
import env
from core.ecan_io import rd_sql, sql_db
import numpy as np
import pandas as pd
from core.classes.hydro import hydro
import rasterio


def get_mean_water_level(min_reading = 2):
    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4) | (well_details_org['WMCRZone'] == 7) |
                                (well_details_org['WMCRZone'] == 8)]  # keep only waimak selwyn and chch zones
    well_details = well_details[well_details['Well_Status'] == 'AE']
    well_details = well_details[pd.notnull(well_details['DEPTH'])]

    screen_details = rd_sql(**sql_db.wells_db.screen_details)
    screen_details.loc[:,'WELL_NO'] = [e.strip() for e in screen_details.loc[:,'WELL_NO']]
    screen_details = screen_details.set_index('WELL_NO')

    well_details.loc[:,'WELL_NO'] = [e.strip() for e in well_details.loc[:,'WELL_NO']]
    well_details = well_details.set_index('WELL_NO')

    data = hydro().get_data(mtypes=['gwl'], sites=list(well_details.index)).data
    data = data.loc['gwl'].reset_index()
    data1 = hydro().get_data(mtypes=['gwl_m'], sites=list(well_details.index)).data
    temp = data1.groupby(level=['mtype', 'site']).describe()[['min', '25%', '50%', '75%', 'mean', 'max', 'count']].round(2)
    sites = list(temp.loc['gwl_m'].index[temp.loc['gwl_m']['count'] >=min_reading])

    data1 = data1.loc['gwl_m',sites].reset_index()
    data = pd.concat((data,data1),axis=0)
    data = data.set_index('site')
    data['month'] = [e.month for e in data.time]
    data['year'] = [e.year for e in data.time]

    rb = rasterio.open(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\LINZ_DEM_8M_Canterbury_2012_31.tif")

    out_data = pd.DataFrame(index=set(data.index))

    for well in out_data.index:
        out_data.loc[well,'nztmx'] = well_details.loc[well,'NZTMX']
        out_data.loc[well,'nztmy'] = well_details.loc[well,'NZTMY']
        out_data.loc[well,'depth'] = well_details.loc[well,'DEPTH']

        ref_level = well_details.loc[well, 'REFERENCE_RL']
        ground_ref_level = well_details.loc[well, 'GROUND_RL']
        if pd.isnull(ground_ref_level):
            ground_ref_level = 0
        if pd.isnull(ref_level):  # if there is no reference level assume it is at the ground from DEM
            x, y = well_details.loc[well, ['NZTMX', 'NZTMY']]
            ref_level = list(rb.sample([(x, y)]))[0][0] + ground_ref_level * -1
            if np.isclose(ref_level, -3.40282306e+38):
                ref_level = np.nan
        out_data.loc[well,'ref_level'] = ref_level
        out_data.loc[well,'ground_level'] = ref_level + ground_ref_level

        if well in screen_details.index:
            top = np.min(screen_details.loc[well,'TOP_SCREEN'])
            bot = np.min(screen_details.loc[well,'BOTTOM_SCREEN'])
            out_data.loc[well,'mid_screen_depth'] = (bot+top)/2
        else:
            out_data.loc[well,'mid_screen_depth'] = well_details.loc[well,'DEPTH'] - 2


    out_data.loc[:,'mid_screen_elv'] = out_data.loc[:,'ground_level'] - out_data.loc[:,'mid_screen_depth']
    data2008 = data.loc[data['year'] >= 2008]
    for val, dat in zip(['_2008',''],[data2008, data]):
        g = dat.loc[(dat.data > -999) & (dat.data < 999)].groupby('site')
        out_data['h2o_dpth{}'.format(val)] = g.aggregate({'data':np.mean})
        out_data['h2o_dpth{}'.format(val)] = g.aggregate({'data':np.min})
        out_data['h2o_dpth{}'.format(val)] = g.aggregate({'data':np.max})
        out_data['h2o_dpth_sd{}'.format(val)] = g.aggregate({'data': np.std})

        temp = dat.loc[(np.in1d(dat['month'],[10,11,12,1,2,3])) & (dat.data > -999) & (dat.data < 999)]
        tempg = temp.groupby('site')
        out_data['h2o_dpth_irr{}'.format(val)] = tempg.aggregate({'data':np.mean})
        out_data['reading_irr{}'.format(val)] = tempg.count().loc[:,'data']

        temp = dat.loc[(np.in1d(dat['month'],[4,5,6,7,8,9])) & (dat.data > -999) & (dat.data < 999)]
        tempg = temp.groupby('site')
        out_data['h2o_dpth_non_irr{}'.format(val)] = tempg.aggregate({'data':np.mean})

        g = dat.groupby('site')
        out_data['readings{}'.format(val)] = g.count().loc[:,'data']

        g = dat.loc[(dat.data <= -999)].groupby('site')
        out_data['count_dry'] = g.count().loc[:,'data']

        g = dat.loc[(dat.data >= 999)].groupby('site')
        out_data['count_flowing'] = g.count().loc[:,'data']

    out_data['h2o_lv_2008'] = out_data.loc[:,'ref_level'] + out_data.loc[:,'h2o_dpth_2008']
    out_data['h2o_lv'] = out_data.loc[:,'ref_level'] + out_data.loc[:,'h2o_dpth']

    out_data['piezo_s_date'] = None
    out_data['seasonal_correction_factor'] = None
    out_data['corrected_wl'] = None
    out_data['corrected_wdpth'] = None
    out_data['rel_weight_fact'] = None

    return out_data


if __name__ == '__main__':
    test = get_mean_water_level()
    test.to_csv(env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/water_levels_for_wells_2_plus_readings.csv"))
    test = get_mean_water_level(1)
    test.to_csv(env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/water_levels_for_wells_1_plus_readings.csv"))



    print 'dpone'