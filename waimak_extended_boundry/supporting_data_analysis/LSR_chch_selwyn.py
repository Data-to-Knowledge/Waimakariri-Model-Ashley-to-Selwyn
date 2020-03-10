# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 1/08/2017 11:36 AM
"""

from __future__ import division
import env
from core.ecan_io import rd_niwa_vcsn
from core.ts.gw.lsrm import poly_import, input_processing, lsrm
from core.misc import save_df
import os
from waimak_extended_boundry import smt
import geopandas as gpd
import pandas as pd
import numpy as np

def run_ds_model(output_results,output_shp):
        #############################################


    #### Parameters

    ### Reading data
    irr_type_dict = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'AQUALINC_NZTM_IRRIGATED_AREA_20160629',
                     'column': 'type'}
    paw_dict = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'LAND_NZTM_NEWZEALANDFUNDAMENTALSOILS',
                'column': 'PAW_MID'}

    bound_shp = env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\model_grid_bound.shp")

    buffer_dis = 10000
    from_date = '2007-08-01'
    to_date = '2016-01-01'
    paw_ratio = 0.67

    ### input data processing
    rain_name = 'rain'
    pet_name = 'pe'

    time_agg = 'W'  # Use 'D' for day, 'W' for week, or 'M' for month
    agg_ts_fun = 'sum'
    grid_res = 200
    crs = 4326
    interp_fun = 'cubic'
    precip_correction = 1.1
    min_irr_area_ratio = 0.01

    irr_mons = [10, 11, 12, 1, 2, 3, 4]

    irr_eff_dict = {'Drip/micro': 1, 'Unknown': 0.8, 'Gun': 0.8, 'Pivot': 0.8, 'K-line/Long lateral': 0.8,
                    'Rotorainer': 0.8, 'Solid set': 0.8, 'Borderdyke': 0.5, 'Linear boom': 0.8, 'Unknown': 0.8,
                    'Lateral': 0.8, 'Wild flooding': 0.5, 'Side Roll': 0.8}
    irr_trig_dict = {'Drip/micro': 0.7, 'Unknown': 0.5, 'Gun': 0.5, 'Pivot': 0.5, 'K-line/Long lateral': 0.5,
                     'Rotorainer': 0.5, 'Solid set': 0.5, 'Borderdyke': 0.5, 'Linear boom': 0.5, 'Unknown': 0.5,
                     'Lateral': 0.5, 'Wild flooding': 0.5, 'Side Roll': 0.5}

    ### Model parameters
    A = 6

    ### Output parameters

    ###########################################
    ### Extract data

    print('Read in the input data')

    irr1, paw1 = poly_import(irr_type_dict, paw_dict, paw_ratio)

    precip_et = rd_niwa_vcsn(mtypes=['precip', 'PET'], sites=bound_shp, buffer_dis=buffer_dis, from_date=from_date,
                             to_date=to_date)

    ##########################################
    ### Process data
    ## Resample met data data

    print('Process the input data for the LSR model')

    model_var, sites_poly = input_processing(precip_et, crs, irr1, paw1, bound_shp, rain_name, pet_name, grid_res,
                                             buffer_dis, interp_fun, agg_ts_fun, time_agg, irr_eff_dict, irr_trig_dict,
                                             min_irr_area_ratio, irr_mons, precip_correction)

    #########################################################
    #### Run the model

    print('Run the LSR model')

    output1 = lsrm(model_var, A)

    #######################################################
    #### Output data

    print('Saving data')

    save_df(output1, output_results, index=False)
    sites_poly.reset_index().to_file(output_shp)



if __name__ == '__main__':
    recalc=False
    out_dir = "{}/m_ex_bd_inputs/DS_LSR".format(smt.sdp)
    output_results = '{}/output_model_domain.csv'.format(out_dir)
    output_shp = '{}/output_model_domain.shp'.format(out_dir)

    if (not os.path.exists(output_results) or not os.path.exists(output_shp)) or recalc:
        run_ds_model(output_results,output_shp)

    shp = gpd.read_file(output_shp)
    results = pd.read_csv(output_results)
    results.loc[:,'time'] = pd.to_datetime(results.loc[:,'time'])
    results = results.loc[(results.time>pd.datetime(2007,12,31))& (results.time<=pd.datetime(2016,12,31))]
    results = results.fillna(0)
    results.loc[:,'rch'] = (results.loc[:,'non_irr_drainage']*(1-results.loc[:,'irr_area_ratio']) +
                          results.loc[:,'irr_drainage']*results.loc[:,'irr_area_ratio'])/results.loc[:,'site_area']/7
    g = results.groupby('site')
    sites = g.aggregate({'rch':np.mean,'paw':np.mean})
    shp2 = pd.merge(shp,sites,left_on='site', right_index=True)
    shp2.to_file("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp), driver='ESRI Shapefile')