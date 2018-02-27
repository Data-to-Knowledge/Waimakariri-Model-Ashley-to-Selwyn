# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 08:03:19 2017

@author: MichaelEK

Functions for the land surface recharge model.
"""
from core.ecan_io.mssql import rd_sql
from core.ts.met.interp import poly_interp_agg
from pandas import merge, concat
from geopandas import read_file
from core.spatial.vector import xy_to_gpd, points_grid_to_poly, spatial_overlays
from numpy import nan, full, arange, seterr, tile, exp


def AET(pet, A, paw_now, paw_max):
    """
    Minhas et al. (1974) function used by David Scott to estimate 'actual ET' from PET and PAW. All inputs must be as floats.
    """

    aet = pet * ((1 - exp(-A*paw_now/paw_max))/(1 - 2*exp(-A) + exp(-A*paw_now/paw_max)))
    return(aet)


def poly_import(irr_type_dict, paw_dict, paw_ratio=0.67):
    """
    Function to import polygon input data. At the moment, these include irrigation type and PAW. Inputs are dictionaries that reference either an MSSQL table with a geometry column or a shapefile. If the dictionary references an sql table then the keys should be 'server', 'database', 'table', and 'column'. If the dictionary references a shapefile, then the keys should be 'shp' and 'column'. All values should be strings.
    """

    if not all([isinstance(irr_type_dict, dict), isinstance(paw_dict, dict)]):
        raise TypeError("'irr_type_dict' and 'paw_dict' must be dictionaries.")

    if 'column' in irr_type_dict.keys():
        if not isinstance(irr_type_dict['column'], str):
            raise TypeError("The key 'column' must be a string.")
    else:
        raise TypeError("The key 'column' must be in the dictionaries.")

    if 'shp' in irr_type_dict.keys():
        if not isinstance(irr_type_dict['shp'], str):
            raise TypeError("If 'shp' is in the dict, then it must be a string path to a shapefile.")
        irr1 = read_file(irr_type_dict['shp'])[[irr_type_dict['column'], 'geometry']]
    elif isinstance(irr_type_dict, dict):
        irr1 = rd_sql(irr_type_dict['server'], irr_type_dict['database'], irr_type_dict['table'], [irr_type_dict['column']], geo_col=True)
    irr1.rename(columns={irr_type_dict['column']: 'irr_type'}, inplace=True)

    if 'shp' in paw_dict.keys():
        if not isinstance(paw_dict['shp'], str):
            raise TypeError("If 'shp' is in the dict, then it must be a string path to a shapefile.")
        paw1 = read_file(paw_dict['shp'])[[paw_dict['column'], 'geometry']]
    elif isinstance(paw_dict, dict):
        paw1 = rd_sql(paw_dict['server'], paw_dict['database'], paw_dict['table'], [paw_dict['column']], geo_col=True)
    paw1.rename(columns={paw_dict['column']: 'paw'}, inplace=True)
    paw1.loc[:, 'paw'] = paw1.loc[:, 'paw'] * paw_ratio

    return(irr1, paw1)


def input_processing(precip_et, crs, irr1, paw1, bound_shp, rain_name, pet_name, grid_res, buffer_dis, interp_fun, agg_ts_fun, time_agg, irr_eff_dict, irr_trig_dict, min_irr_area_ratio=0.01, irr_mons=[10, 11, 12, 1, 2, 3, 4], precip_correction=1.1):
    """
    Function to process the input data for the lsrm. Outputs a DataFrame of the variables for the lsrm.
    """
    seterr(invalid='ignore')

    ## Load and resample precip and et
    bound = read_file(bound_shp)

    new_rain = poly_interp_agg(precip_et, crs, bound_shp, rain_name, 'time', 'x', 'y', buffer_dis, grid_res, grid_res, interp_fun=interp_fun, agg_ts_fun=agg_ts_fun, period=time_agg) * precip_correction
    new_rain.name = 'precip'

    new_et = poly_interp_agg(precip_et, crs, bound_shp, pet_name, 'time', 'x', 'y', buffer_dis, grid_res, grid_res, interp_fun=interp_fun, agg_ts_fun=agg_ts_fun, period=time_agg)
    new_et.name = 'pet'

    new_rain_et = concat([new_rain, new_et], axis=1)

    ## convert new point locations to geopandas
    time1 = new_rain_et.index.levels[0][0]
    grid1 = new_rain_et.loc[time1].reset_index()[['x', 'y']]
    grid2 = xy_to_gpd(grid1.index, 'x', 'y', grid1, bound.crs)
    grid2.columns = ['site', 'geometry']

    all_times = new_rain_et.index.levels[0]
    new_rain_et.loc[:, 'site'] = tile(grid1.index, len(all_times))

    ## Convert points to polygons
    sites_poly = points_grid_to_poly(grid2, 'site')

    ## process polygon data
    # Select polgons within boundary

    sites_poly_union = sites_poly.unary_union
    irr2 = irr1[irr1.intersects(sites_poly_union)]
    irr3 = irr2[irr2.irr_type.notnull()]
    paw2 = paw1[paw1.intersects(sites_poly_union)]
    paw3 = paw2[paw2.paw.notnull()]

    # Overlay intersection
    sites_poly1 = spatial_overlays(sites_poly, bound, how='intersection')[['site', 'geometry']]
    sites_poly2 = sites_poly1.dissolve('site')
    sites_poly2.crs = sites_poly.crs
    sites_poly_area = sites_poly2.area.round(2)
    sites_poly3 = sites_poly2.reset_index()

    irr4 = spatial_overlays(irr3, sites_poly3, how='intersection')
    paw4 = spatial_overlays(paw3, sites_poly3, how='intersection')

    irr4['area'] = irr4.geometry.area.round()
    irr5 = irr4[irr4.area >= 1].drop(['idx1', 'idx2'], axis=1).copy()

    paw4['area'] = paw4.geometry.area.round()
    paw5 = paw4.loc[(paw4.area >= 1)].drop(['idx1', 'idx2'], axis=1).copy()
    paw5.loc[paw5.paw <= 0, 'paw'] = 1

    # Add in missing PAW values - Change later to something more useful if needed
    mis_sites_index = ~sites_poly3.site.isin(paw5.site)
    sites_poly3['area'] = sites_poly3.area.round()

    paw6 = concat([paw5, sites_poly3[mis_sites_index]])
    paw6.loc[paw6.paw.isnull(), 'paw'] = 1

    # Aggregate by site weighted by area to estimate a volume
    paw_area1 = paw6[['paw', 'site', 'area']].copy()
    paw_area1.loc[:, 'paw_vol'] = paw_area1['paw'] * paw_area1['area']
    paw7 = ((paw_area1.groupby('site')['paw_vol'].sum() / paw_area1.groupby('site')['area'].sum()) * sites_poly_area * 0.001).round(2)

    site_irr_area = irr5.groupby('site')['area'].sum()
    irr_eff1 = irr5.replace({'irr_type': irr_eff_dict})
    irr_eff1.loc[:, 'irr_eff'] = irr_eff1['irr_type'] * irr_eff1['area']
    irr_eff2 = (irr_eff1.groupby('site')['irr_eff'].sum() / site_irr_area).round(3)

    irr_trig1 = irr5.replace({'irr_type': irr_trig_dict})
    irr_trig1.loc[:, 'irr_trig'] = irr_trig1['irr_type'] * irr_trig1['area']
    irr_trig2 = (irr_trig1.groupby('site')['irr_trig'].sum() / site_irr_area).round(3)

    irr_area_ratio1 = (site_irr_area/sites_poly_area).round(3)

    poly_data1 = concat([paw7, sites_poly_area, irr_eff2, irr_trig2, irr_area_ratio1], axis=1)
    poly_data1.columns = ['paw', 'site_area', 'irr_eff', 'irr_trig', 'irr_area_ratio']
    poly_data1.loc[poly_data1['irr_area_ratio'] < min_irr_area_ratio, ['irr_eff', 'irr_trig', 'irr_area_ratio']] = nan

    ## Combine time series with polygon data
    new_rain_et1 = new_rain_et[new_rain_et['site'].isin(sites_poly2.index)]

    input1 = merge(new_rain_et1.reset_index(), poly_data1.reset_index(), on='site', how='left')

    ## Convert precip and et to volumes
    input1.loc[:, ['precip', 'pet']] = (input1.loc[:, ['precip', 'pet']].mul(input1.loc[:, 'site_area'], axis=0) * 0.001).round(2)

    ## Remove irrigation parameters during non-irrigation times
    input1.loc[~input1.time.dt.month.isin(irr_mons), ['irr_eff', 'irr_trig']] = nan

    ## Run checks on the input data

#    print('Running checks on the prepared input data')

    null_time = input1.loc[input1.time.isnull(), 'time']
    null_x = input1.loc[input1.x.isnull(), 'x']
    null_y = input1.loc[input1.y.isnull(), 'y']
    null_pet = input1.loc[input1['pet'].isnull(), 'pet']
    null_rain = input1.loc[input1['precip'].isnull(), 'precip']
    null_paw = input1.loc[input1.paw.isnull(), 'paw']
    not_null_irr_eff = input1.loc[input1.irr_eff.notnull(), 'irr_eff']

    if not null_time.empty:
        raise ValueError('Null values in the time variable')
    if not null_x.empty:
        raise ValueError('Null values in the x variable')
    if not null_y.empty:
        raise ValueError('Null values in the y variable')
    if not null_pet.empty:
        raise ValueError('Null values in the pet variable')
    if not null_rain.empty:
        raise ValueError('Null values in the rain variable')
    if not null_paw.empty:
        raise ValueError('Null values in the paw variable')
    if not_null_irr_eff.empty:
        raise ValueError('No values for irrigation variables')

    if input1['time'].dtype.name != 'datetime64[ns]':
        raise ValueError('time variable must be a datetime64[ns] dtype')
    if input1['x'].dtype != float:
        raise ValueError('x variable must be a float dtype')
    if input1['y'].dtype != float:
        raise ValueError('y variable must be a float dtype')
    if input1['pet'].dtype != float:
        raise ValueError('pet variable must be a float dtype')
    if input1['precip'].dtype != float:
        raise ValueError('precip variable must be a float dtype')
    if input1['paw'].dtype != float:
        raise ValueError('paw variable must be a float dtype')
    if input1['irr_eff'].dtype != float:
        raise ValueError('irr_eff variable must be a float dtype')
    if input1['irr_trig'].dtype != float:
        raise ValueError('irr_trig variable must be a float dtype')
    if input1['irr_area_ratio'].dtype != float:
        raise ValueError('irr_area_ratio variable must be a float dtype')

    ## Return dict
    return(input1, sites_poly2)


def lsrm(model_var, A, include_irr=True):
    """
    The lsrm.
    """
    seterr(invalid='ignore')

    ## Make initial variables
    all_times = model_var['time'].unique()
    sites = model_var['site'].unique()

    ## Prepare individual variables
    # Input variables
    paw_val = model_var['paw'].values
    irr_area_ratio_val = model_var['irr_area_ratio'].copy()
    irr_area_ratio_val[irr_area_ratio_val.isnull()] = 0

    irr_paw_val = paw_val * irr_area_ratio_val.values
    non_irr_paw_val = paw_val - irr_paw_val
    irr_paw_val[irr_paw_val <= 0 ] = nan
    non_irr_paw_val[non_irr_paw_val <= 0 ] = nan

    irr_rain_val = model_var['precip'].values * irr_area_ratio_val.values
    non_irr_rain_val = model_var['precip'].values - irr_rain_val

    irr_pet_val = model_var['pet'].values * irr_area_ratio_val.values
    non_irr_pet_val = model_var['pet'].values - irr_pet_val

    irr_eff_val = model_var['irr_eff'].values
    irr_trig_val = model_var['irr_trig'].values

    time_index = arange(len(model_var)).reshape(len(all_times), len(sites))

    # Initial conditions
    w_irr = irr_paw_val[time_index[0]].copy()
    w_non_irr = non_irr_paw_val[time_index[0]].copy()

    # Output variables
    out_w_irr = full(len(irr_eff_val), nan)
    out_irr_demand = out_w_irr.copy()
    out_w_non_irr = out_w_irr.copy()
    out_irr_aet = out_w_irr.copy()
    out_non_irr_aet = out_w_irr.copy()

    out_irr_drainage = out_w_irr.copy()
    out_non_irr_drainage = out_w_irr.copy()

    ## Run the model
    for i in time_index:

        if include_irr:
            ### Irrigation bucket
            i_irr_paw = irr_paw_val[i]

            ## Calc AET and perform the initial water balance
            irr_aet_val = AET(irr_pet_val[i], A, w_irr, i_irr_paw)
            out_irr_aet[i] = irr_aet_val.copy()
            w_irr = w_irr + irr_rain_val[i] - irr_aet_val

            ## Check and calc the GW drainage from excessive precip
            irr_drainge_bool = w_irr > i_irr_paw
            if any(irr_drainge_bool):
                temp_irr_draiange = w_irr[irr_drainge_bool] - i_irr_paw[irr_drainge_bool]
                out_irr_drainage[i[irr_drainge_bool]] = temp_irr_draiange
                w_irr[irr_drainge_bool] = i_irr_paw[irr_drainge_bool]
            out_w_irr[i] = w_irr.copy()

            ## Check and calc drainage from irrigation if w is below trigger
            irr_paw_ratio = w_irr/i_irr_paw
            irr_trig_bool = irr_paw_ratio <= irr_trig_val[i]
            if any(irr_trig_bool):
                diff_paw = i_irr_paw[irr_trig_bool] - w_irr[irr_trig_bool]
                out_irr_demand[i[irr_trig_bool]] = diff_paw.copy()
                irr_drainage = diff_paw/irr_eff_val[i][irr_trig_bool] - diff_paw
                out_irr_drainage[i[irr_trig_bool]] = irr_drainage.copy()
                w_irr[irr_trig_bool] = i_irr_paw[irr_trig_bool].copy()

            out_w_irr[i] = w_irr.copy()

        ### Non-irrigation bucket
        i_non_irr_paw = non_irr_paw_val[i]

        ## Calc AET and perform the initial water balance
        non_irr_aet_val = AET(non_irr_pet_val[i], A, w_non_irr, i_non_irr_paw)
        out_non_irr_aet[i] = non_irr_aet_val.copy()
        w_non_irr = w_non_irr + non_irr_rain_val[i] - non_irr_aet_val
        w_non_irr[w_non_irr < 0] = 0

        ## Check and calc the GW drainage from excessive precip
        non_irr_drainge_bool = w_non_irr > i_non_irr_paw
        if any(non_irr_drainge_bool):
            temp_non_irr_draiange = w_non_irr[non_irr_drainge_bool] - i_non_irr_paw[non_irr_drainge_bool]
            out_non_irr_drainage[i[non_irr_drainge_bool]] = temp_non_irr_draiange.copy()
            w_non_irr[non_irr_drainge_bool] = i_non_irr_paw[non_irr_drainge_bool]
        out_w_non_irr[i] = w_non_irr.copy()

    ### Put results into dataframe

    output1 = model_var.copy()
    output1.loc[:, 'non_irr_aet'] = out_non_irr_aet.round(2)
    if include_irr:
        output1.loc[:, 'irr_aet'] = out_irr_aet.round(2)
        output1.loc[:, 'irr_paw'] = irr_paw_val.round(2)
        output1.loc[:, 'w_irr'] = out_w_irr.round(2)
        output1.loc[:, 'irr_drainage'] = out_irr_drainage.round(2)
        output1.loc[:, 'irr_demand'] = out_irr_demand.round(2)
    else:
        output1.drop(['irr_eff', 'irr_trig', 'irr_area_ratio'], axis=1, inplace=True)
    output1.loc[:, 'non_irr_paw'] = non_irr_paw_val.round(2)
    output1.loc[:, 'w_non_irr'] = out_w_non_irr.round(2)
    output1.loc[:, 'non_irr_drainage'] = out_non_irr_drainage.round(2)

    ### Return output dataframe
    return(output1)


























