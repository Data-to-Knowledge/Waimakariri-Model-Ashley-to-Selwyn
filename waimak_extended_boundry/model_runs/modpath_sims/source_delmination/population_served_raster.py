# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 9:47 AM

script to create a map of the population served by a given area
"""

from __future__ import division
import env
import os
from waimak_extended_boundry import \
    get_well_ids
from glob import glob
from waimak_extended_boundry import smt
import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter


# below from \\gisdata\ProjectArchive\SCI\2015_2016\EMG\2015_2016\Waimakariri\Groundwater\Groundwater Quantity\Consented GW abstraction\WDC Flow & Population Data.xlsx
wdc_population = {'Cust':350,
 'Fernside':85*3, # do not have an estimate assumed 3 people on average for all of the connections
 'Kaiapoi':10000,
 'Kairaki':200, # do not have an estimate, likely a back up well, we just threw in 200 ppl
 'Mandeville':1800,
 'Ohoka':233,
 'Oxford Urban':3315, # sum of oxford urban and rural
 'Pegasus':1050,
 'Poyntzs Road':200,
 'Rangiora':13500,
 'Waikuku':1150,
 'West Eyreton':120
                  }


def calc_pop_served_raster(num_private, n_threshold=0):
    """
    creates a raster with the population served
    :param num_private: the assumed number of people served by each private well
    :param n_threshold: the number below which we do not add the population to the raster default 0 all wells added
    :return:
    """
    base_shp_dir = env.sci("Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results"
                           "\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons")
    outdata = smt.get_empty_model_grid()
    # wdc
    wdc_n = pd.read_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\wdc_wells\stocastic_n_wdc_wells\stocastic_set_with_trans_overview.csv",
        index_col=0)
    wdc_shp_dir = os.path.join(base_shp_dir, "wdc_wells")
    wdc_shps = glob(os.path.join(wdc_shp_dir, '*.shp'))
    for shp in wdc_shps:
        nm = os.path.basename(shp).split('.')[0]
        if wdc_n.loc['wdc_{}'.format(nm), '50%_cmp_load'] < n_threshold:
            continue
        temp = np.isfinite(smt.shape_file_to_model_array(shp, 'Id', True)).astype(int)
        outdata += temp * wdc_population[nm]

    private_n = pd.read_csv(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\private_wells\stocastic_n_private_wells\stocastic_set_with_trans_overview.csv",
        index_col=0)
    # private wells
    private_shp_dir = os.path.join(base_shp_dir, "private_wells")
    well_ids = get_well_ids()
    private_shps = glob(os.path.join(private_shp_dir, '*.shp'))
    for shp in private_shps:
        nm = os.path.basename(shp).split('.')[0]
        if private_n.loc['{}'.format(nm), '50%_cmp_load'] < n_threshold:
            continue
        temp = np.isfinite(smt.shape_file_to_model_array(shp, 'Id', True)).astype(int)
        num_ser = len(well_ids.loc[(well_ids.Zone_1 == nm) | (well_ids.zone_2 == nm)])
        outdata += temp * num_ser * num_private
    return outdata


if __name__ == '__main__':
    blur_factor = .5
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\population_served"
    test = calc_pop_served_raster(3).astype(float)
    test[np.isclose(test, 0)] = np.nan
    test = gaussian_filter(test, blur_factor)
    smt.array_to_raster(os.path.join(outdir, 'population_served_blurred.tif'), test, 0)

    for n_limit in [3.5, 5.6, 11.3]:
        test2 = calc_pop_served_raster(3, n_threshold=n_limit).astype(float)
        test2[np.isclose(test2, 0)] = np.nan
        test2 = gaussian_filter(test2, blur_factor)
        smt.array_to_raster(os.path.join(outdir, 'population_served_n_over_{}_blurred.tif'.format(n_limit)), test2, 0)
