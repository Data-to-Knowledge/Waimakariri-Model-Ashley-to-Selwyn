# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/09/2017 1:56 PM
"""

from __future__ import division
import env
import pandas as pd
import numpy as np
from glob import glob
import os

test_path = r"K:\niwa_netcdf\lsrm\lsrm_results\vcsn_65perc.h5"


def gen_water_year_average_lsr_irr(path):
    org_data = pd.read_hdf(path).drop(['x', 'y','irr_eff', 'irr_trig','non_irr_aet', 'irr_aet','irr_paw', 'w_irr'], axis=1, errors='ignore')
    if 'irr_drainage' not in org_data.keys():
        org_data.loc[:, 'irr_drainage'] = 0
    if 'irr_demand' not in org_data.keys():
        org_data.loc[:, 'irr_demand'] = 0
    if 'irr_area_ratio' not in org_data.keys():
        org_data.loc[:, 'irr_area_ratio'] = 0

    org_data.loc[org_data.irr_drainage.isnull(), 'irr_drainage'] = 0
    org_data.loc[org_data.non_irr_drainage.isnull(), 'non_irr_drainage'] = 0
    org_data['total_drainage'] = org_data.irr_drainage + org_data.non_irr_drainage
    g = org_data.groupby([pd.Grouper('site'), pd.Grouper(key='time', freq='A-JUN')])
    temp = org_data.drop_duplicates('site').set_index('site')
    temp.loc[temp.irr_area_ratio.isnull(), 'irr_area_ratio'] = 0
    outdata = g.aggregate(
        {'precip': np.nansum, 'pet': np.nansum, 'paw': np.nansum, 'irr_drainage': np.nansum, 'irr_demand': np.nansum,
         'non_irr_drainage': np.nansum, 'total_drainage': np.nansum}).reset_index()
    outdata['site_area'] = outdata.loc[:, 'site']
    outdata['irr_area_ratio'] = outdata.loc[:, 'site']
    outdata = outdata.replace({'site_area': temp.site_area.to_dict()})
    outdata = outdata.replace({'irr_area_ratio': temp.irr_area_ratio.to_dict()})
    for key in ['precip', 'pet', 'paw', 'total_drainage']:  # raising memory error so trying to do it one at a time
        outdata[key] *= (1 / (365 * outdata.site_area))
    outdata['irr_drainage'] *= (1 / (365 * outdata.site_area * outdata.irr_area_ratio))  # per area of irrigation
    outdata['irr_demand'] *= (1 / (365 * outdata.site_area * outdata.irr_area_ratio))  # m3/day
    outdata['non_irr_drainage'] = (
    1 / (365 * outdata.site_area * (1 - outdata.irr_area_ratio)))  # per area of non-irrigation
    return outdata


def gen_all_wy_averages(base_dir):
    """
    generate all the water year mean recharge needs to be run with 20+ g of ram
    :param base_dir: the directory that holds all of the LSRM results
    :return:
    """
    out_dir = '{}/water_year_means'.format(base_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    paths = glob('{}/*.h5'.format(base_dir))
    for path in paths:
        print (path)
        outdata = gen_water_year_average_lsr_irr(path)
        outdata.to_hdf(os.path.join(out_dir, 'wym_{}'.format(os.path.basename(path))), 'wym', mode='w')


if __name__ == '__main__':
    gen_all_wy_averages(env.gw_met_data("niwa_netcdf/lsrm/lsrm_results"))
