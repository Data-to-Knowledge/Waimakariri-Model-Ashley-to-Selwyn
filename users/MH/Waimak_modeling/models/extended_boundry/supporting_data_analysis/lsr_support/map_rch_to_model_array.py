# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/09/2017 5:13 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from warnings import warn
from copy import deepcopy


def map_rch_to_array(hdf, method, period_center, mapping_shp, period_length=10, return_irr_demand=False,
                     rch_quanity='total_drainage',site_list=None):
    """
    takes a hdf and maps it into an array for the model, assumes that the hdf file is a water year amalgimation
    :param hdf: the hdf file that contains the water year average of the LSR from the LSRM
    :param method: one of:
                            'mean': mean of the full dataset
                            'period_mean': mean of the period defined above if the period is longer than the data set
                                           this is equivelent to 'mean' but warning is passed
                            '3_lowest_con_mean': mean of the three lowest consecuative years for sites defined by site_list
                            'lowest_year': simply the lowest year for sites defined by site_list
    :param mapping_shp: the path to the shapefile that links the data to the geometry
    :param period_center: the year of the center (int)
    :param period_length: the length of the period in years (int) must be even
    :param return_irr_demand: if true then returen both recharge and average annual (spatially summed) irrigation demand (single number)
    :param rch_quanity: the value to use
    :param site_list: a list of site number (from the data) which should be considered for lowest year.  if None then all considered
    :return: array of the shape of (smt.rows, smt.cols) if return irr_demand then (rch, irr_demand)
    """
    # quick check of inputs
    if method not in ['mean', 'period_mean', '3_lowest_con_mean', 'lowest_year']:
        raise ValueError('unsuported method')

    if not isinstance(period_center, int) and period_center is not None:
        raise TypeError('period_center must be int or None')

    if not isinstance(period_length, int) and period_length is not None:
        raise TypeError('period_length must be int or None')

    if isinstance(period_length, int):
        if period_length % 2 != 0:
            raise ValueError('period length must be even')

    outdata = smt.shape_file_to_model_array(mapping_shp, 'site', True)
    data = pd.read_hdf(hdf).reset_index()
    data.loc[:, 'year'] = [e.year for e in data.time]
    if method == 'mean':
        map_data = data.groupby('site').aggregate({rch_quanity: np.mean})
        if return_irr_demand:
            ird_map_data = data.groupby('site').aggregate({'irr_demand': np.mean})
    else:
        data = data.loc[
            (data.year >= period_center - period_length / 2) & (data.year <= period_center + period_length / 2)]

        if site_list is None:
            year_data = data.groupby('year').aggregate({rch_quanity: np.sum, 'irr_demand': np.sum})
        else:
            if not set(site_list).issubset(set(data.site)):
                raise ValueError('sites not defined: {}'.format(set(site_list)-set(data.site)))
            year_data = data.loc[np.in1d(data.site,site_list)].groupby('year').aggregate({rch_quanity: np.sum, 'irr_demand': np.sum})

        expected_years = range(int(period_center - period_length / 2), int(period_center + period_length / 2 + 1))
        if not all(np.in1d(expected_years, year_data.index)):
            warn('not all period years present taking mean of full dataset, years missing: {}'.format(
                set(expected_years) - set(year_data.index)))

        if method == 'period_mean':
            map_data = data.groupby('site').aggregate({rch_quanity: np.mean})
            if return_irr_demand:
                ird_map_data = data.groupby('site').aggregate({'irr_demand': np.mean})
        elif method == 'lowest_year':
            low_year = year_data[rch_quanity].argmin()
            map_data = data.loc[data.year == low_year].groupby('site').aggregate({rch_quanity: np.mean})
            if return_irr_demand:
                ird_map_data = data.loc[data.year == low_year].groupby('site').aggregate({'irr_demand': np.mean})
        elif method == '3_lowest_con_mean':
            if len(year_data) < 3:
                raise ValueError('3 year consectuative mean cannot be calculated on data with less than 3 years')
            elif len(year_data) == 3:
                warn('only 3 years of data, same as period mean')
            low_year = year_data.rolling(3).sum()[rch_quanity].argmin()
            low_years = [low_year - e for e in [2, 1, 0]]
            map_data = data.loc[np.in1d(data.year, low_years)].groupby('site').aggregate({rch_quanity: np.mean})
            if return_irr_demand:
                ird_map_data = data.loc[np.in1d(data.year, low_years)].groupby('site').aggregate({'irr_demand': np.mean})
        else:
            raise ValueError('should not get here')
    outdata[np.isnan(outdata)] = -999
    outdata_ird = deepcopy(outdata)
    temp_dict = map_data[rch_quanity].to_dict()
    temp_dict[-999] = np.nan
    outdata = vec_translate(outdata, temp_dict)
    if return_irr_demand:
        ird_temp_dict = ird_map_data['irr_demand'].to_dict()
        ird_temp_dict[-999] = np.nan
        ird_outdata = vec_translate(outdata_ird, ird_temp_dict)

    if return_irr_demand:
        return outdata, ird_outdata
    else:
        return outdata


# quickest way to assign is to loop or to use below
def vec_translate(a, d):
    return np.vectorize(d.__getitem__)(a)


if __name__ == '__main__':
    from glob import glob
    import matplotlib.pyplot as plt

    test_path = r"K:\niwa_netcdf\lsrm\lsrm_results\vcsn_65perc.h5"
    hdf_paths = glob(r"K:\niwa_netcdf\lsrm\lsrm_results\water_year_means\wym_vcsn_*perc.h5")
    map_shp = r"K:\niwa_netcdf\lsrm\lsrm_results\test\output_test2.shp"
    methods = ['mean', 'period_mean', '3_lowest_con_mean', 'lowest_year']

    rch_arrays = {}
    ird_demand = {}
    import time
    t = time.time()
    for hdf_path in hdf_paths:
        print (hdf_path)
        per = int(hdf_path.split('_')[-1].split('.')[0].strip('perc'))
        rch_arrays[per], ird_demand[per] = map_rch_to_array(hdf=hdf_path,
                                                            method=methods[0],
                                                            period_center=2012,
                                                            mapping_shp=map_shp,
                                                            period_length=10,
                                                            return_irr_demand=True)
    print(time.time() - t)
    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = np.nan
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    w_idx = np.isclose(zones, 4)

    print('lsr dif 50 - 100 waimak')
    print (np.nansum((rch_arrays[50] - rch_arrays[100])[w_idx] * 200 * 200) / 86400)
    print ('lsr dif 50 - 100 all')
    print (np.nansum((rch_arrays[50] - rch_arrays[100]) * 200 * 200) / 86400)
    print ('lsr dif 65 - 100 waimak')
    print (np.nansum((rch_arrays[65] - rch_arrays[100])[w_idx] * 200 * 200) / 86400)
    print ('lsr dif 65 - 100 all')
    print (np.nansum((rch_arrays[65] - rch_arrays[100]) * 200 * 200) / 86400)
    print ('lsr dif 80 - 100 waimak')
    print (np.nansum((rch_arrays[80] - rch_arrays[100])[w_idx] * 200 * 200) / 86400)
    print ('lsr dif 80 - 100 all')
    print (np.nansum((rch_arrays[80] - rch_arrays[100]) * 200 * 200) / 86400)
    smt.plt_matrix(rch_arrays[50] / rch_arrays[100], title='50/100')
    smt.plt_matrix(rch_arrays[65] / rch_arrays[100], title='65/100')
    smt.plt_matrix(rch_arrays[80] / rch_arrays[100], title='80/100')
    smt.plt_matrix(rch_arrays[80], title='80')

    plt.show()