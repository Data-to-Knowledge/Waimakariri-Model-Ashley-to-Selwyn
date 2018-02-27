# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/10/2017 1:33 PM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import numpy as np
import pandas as pd
from glob import glob
import os
import geopandas as gpd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import itertools

zone_lats = {}
zone_lons = {}

bounds = gpd.read_file(
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\lsr_test_zones_wgs.shp")
for _id, geom, zone, in bounds.itertuples(False, None):
    zone_lats[zone] = geom.bounds[1], geom.bounds[3]
    zone_lons[zone] = geom.bounds[0], geom.bounds[2]


def make_save_rcppast_year_amalg(variable, amalg_type, outpath, groupby=('year')):
    """

    :param variable:
    :param amalg_type: e.g np.sum
    :return:
    """
    if variable == 'pe':
        paths = glob(env.gw_met_data(
            r"niwa_netcdf\climate_projections\RCPpast\*\PE_VCSN_*_RCPpast_1971_2005_south-island_p05_daily_ECan.nc"))
        at = np.nanmean
    elif variable == 'rain':
        paths = glob(env.gw_met_data(
            r"niwa_netcdf\climate_projections\RCPpast\*\TotalPrecipCorr_VCSN_*_RCPpast_1971_2005_south-island_p05_daily_ECan.nc"))
        at = np.nansum
    outdata = []
    paths.append(r"Y:\VirtualClimate\vcsn_precip_et_2016-06-06.nc")
    for path in paths:
        model = os.path.basename(os.path.dirname(path))
        data = nc.Dataset(path)
        lat, lon = np.array(data['latitude']), np.array(data['longitude'])
        for zone in zone_lats.keys():
            lats = zone_lats[zone]
            lons = zone_lons[zone]
            lat_idx = np.where((lat >= lats[0]) & (lat <= lats[1]))
            lon_idx = np.where((lon >= lons[0]) & (lon <= lons[1]))
            temp_var = np.array(data[variable])
            temp_var[np.isclose(temp_var, data.variables[variable]._FillValue)] = np.nan
            print(np.isnan(temp_var).sum() / len(temp_var.flatten()))
            if path == r"Y:\VirtualClimate\vcsn_precip_et_2016-06-06.nc":
                temp_var = amalg_type(temp_var[lon_idx[0], :, :][:, lat_idx[0], :],
                                      axis=(1, 0))
            else:
                temp_var = amalg_type(temp_var[:, lat_idx[0]][:, :, lon_idx[0]],
                                      axis=(2, 1))
            years = [e.year for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
            month = [e.month for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
            temp_outdata = pd.DataFrame({'year': years, 'month': month, variable: temp_var})
            temp_outdata.loc[:, 'model'] = model
            temp_outdata.loc[:, 'zone'] = zone
            temp_outdata = temp_outdata.groupby(list(groupby) + ['model', 'zone']).aggregate(
                {variable: at}).reset_index()
            outdata.append(temp_outdata)

    outdata = pd.concat(outdata)
    outdata.to_csv(outpath)



def plt_data_ts(plt_idx=['year'], title=None):
    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore"

    if plt_idx == ['year', 'month']:
        path = os.path.join(base_dir, 'year_month_zonal_pe_comp.csv')
        data = pd.read_csv(path)
        data.loc[:, 'datetime'] = [pd.datetime(y, m, 15) for y, m in
                                   data.loc[:, ['year', 'month']].itertuples(False, None)]
        x = 'datetime'
    elif plt_idx == ['year']:
        path = os.path.join(base_dir, 'year_zonal_pe_comp.csv')
        data = pd.read_csv(path)
        x = 'year'
    elif plt_idx == ['month']:
        path = os.path.join(base_dir, 'month_zonal_pe_comp.csv')
        data = pd.read_csv(path)
        x = 'month'
    else:
        raise ValueError('plt idx not implemented')

    g = sns.FacetGrid(data, row='zone',legend_out=False)
    g.map_dataframe(sns.pointplot, x, 'pe', 'model', palette=sns.color_palette("Set1", n_colors=8, desat=.5))
    g.add_legend()
    g.fig.suptitle(title)
    return g


def make_data():
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore"
    outpath = os.path.join(outdir, 'year_zonal_pe_comp.csv')
    make_save_rcppast_year_amalg('pe', np.nanmean, outpath, ['year'])

    outpath = os.path.join(outdir, 'month_zonal_pe_comp.csv')
    make_save_rcppast_year_amalg('pe', np.nanmean, outpath, ['month'])

    outpath = os.path.join(outdir, 'year_month_zonal_pe_comp.csv')
    make_save_rcppast_year_amalg('pe', np.nanmean, outpath, ['year', 'month'])

def plot_all_data_ts():
    g1 = plt_data_ts(['year'], 'year')
    g2 = plt_data_ts(['year', 'month'], 'month-year')
    g3 = plt_data_ts(['month'], 'month')
    plt.show()

def plot_month_year_relationship(outpath):
    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore"
    path = os.path.join(base_dir, 'year_month_zonal_pe_comp.csv')
    data = pd.read_csv(path)
    data.loc[:, 'datetime'] = [pd.datetime(y, m, 15) for y, m in
                               data.loc[:, ['year', 'month']].itertuples(False, None)]
    vcsn_data = data.loc[data.model=='VirtualClimate']
    rcm_data = data.loc[~(data.model=='VirtualClimate')]
    plot_data = pd.merge(rcm_data,vcsn_data,on=['zone','datetime'])
    # do regression
    outdata = pd.DataFrame(data=list(itertools.product(set(plot_data.zone),set(plot_data.model_x))), columns=['zone','model'])
    for key in ['slope', 'intercept', 'r_value', 'p_value', 'std_err']:
        outdata[key] = np.nan
    plot_data.loc[:,'pe_x'] = np.log(plot_data.pe_x) # the end looks like log it better
    plot_data.loc[:,'pe_y'] = np.log(plot_data.pe_y)
    title = 'log_transform'

    for i,zone, model_x in outdata.loc[:,['zone','model']].itertuples(True,None):
        temp = plot_data.loc[(plot_data.model_x==model_x) & (plot_data.zone == zone)]
        slope, intercept, r_value, p_value, std_err = stats.linregress(temp['pe_x'], temp['pe_y'])
        outdata.loc[i, ['slope', 'intercept', 'r_value', 'p_value', 'std_err']] = slope, intercept, r_value, p_value, std_err
    outdata.to_csv(outpath)
    # plot
    g = sns.FacetGrid(plot_data, row='zone', col='model_x')
    g.map_dataframe(sns.regplot, x='pe_x', y='pe_y')
    g.fig.suptitle(title)
    g2 = sns.FacetGrid(plot_data, row='zone', col='model_x')
    g2.map_dataframe(sns.residplot, x='pe_x', y='pe_y')
    g2.fig.suptitle(title)
    plt.show()




if __name__ == '__main__':
    #make_data()
    #plot_all_data_ts()
    plot_month_year_relationship(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\test_regressions.csv")
