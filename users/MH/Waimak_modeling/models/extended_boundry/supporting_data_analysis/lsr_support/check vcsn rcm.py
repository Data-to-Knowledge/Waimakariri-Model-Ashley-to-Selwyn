# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 10/10/2017 5:47 PM
"""

from __future__ import division
from core import env
import netCDF4 as nc
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns
import os

lats = (-43.474998, -43.275002)
lons = (171.975006, 172.725006)


def make_save_rcppast_year_amalg(variable, amalg_type, groupby=['year'],outpath=None):
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
    for path in paths:
        # heaps of these data are nan why? because they are in the ocean
        data = nc.Dataset(path)
        lat, lon = np.array(data['latitude']), np.array(data['longitude'])
        lat_idx = np.where((lat >= lats[0]) & (lat <= lats[1]))
        lon_idx = np.where((lon >= lons[0]) & (lon <= lons[1]))
        temp_var = np.array(data[variable])
        temp_var[np.isclose(temp_var, data.variables[variable]._FillValue)] = np.nan
        print(np.isnan(temp_var).sum() / len(temp_var.flatten()))
        temp_var = amalg_type(temp_var[:, lat_idx[0]][:, :, lon_idx[0]],
                              axis=(2, 1))
        outdata.append(temp_var[:, np.newaxis])
    outdata = np.concatenate(outdata, axis=1).mean(axis=1)
    years = [e.year for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
    month = [e.month for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
    outdata = pd.DataFrame({'year': years, 'month': month, variable: outdata})
    outdata = outdata.groupby(groupby).aggregate({variable: at})
    print('done')
    if outpath is not None:
        outdata.to_csv(outpath)
    return outdata


def make_save_vcsn_year_mean(variable, amalg_type, groupby=['year'], outpath=None):
    data_path = "C:\Users\MattH\Desktop\{}_rcm_year.csv".format(variable)
    # units, pe is in mm/day and rain is in mm
    if variable == 'pe':
        at = np.nanmean
    elif variable == 'rain':
        at = np.nansum
    path = r"\\fileservices02\managedshares\data\VirtualClimate\vcsn_precip_et_2016-06-06.nc"
    data = nc.Dataset(path)
    lat, lon = np.array(data['latitude']), np.array(data['longitude'])
    lat_idx = np.where((lat >= lats[0]) & (lat <= lats[1]))
    lon_idx = np.where((lon >= lons[0]) & (lon <= lons[1]))
    temp_var = np.array(data[variable])
    temp_var[np.isclose(temp_var, data.variables[variable]._FillValue)] = np.nan
    print(np.isnan(temp_var).sum() / len(temp_var.flatten()))
    temp_var = amalg_type(temp_var[lon_idx[0], :, :][:, lat_idx[0], :],
                          axis=(0, 1))
    years = [e.year for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
    month = [e.month for e in nc.num2date(np.array(data.variables['time']), data.variables['time'].units)]
    outdata = pd.DataFrame({'year': years, 'month': month, variable: temp_var})
    outdata = outdata.groupby(groupby).aggregate({variable: at})
    if outpath is not None:
        outdata.to_csv(outpath)
    return outdata


if __name__ == '__main__':
    # rainfig, rainax = plt.subplots()
    # rcp_past_rain = make_save_rcppast_year_amalg('rain', np.nanmean)
    # vcsn_rain = make_save_vcsn_year_mean('rain',np.nanmean)
    # rainax.plot(rcp_past_rain.index,rcp_past_rain.rain,'r',label='rcpPAST rain')
    # rainax.plot(vcsn_rain.index,vcsn_rain.rain,'b',label='vcsn rain')
    # rainax.legend()
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore"
    pefig, peax = plt.subplots()
    vcsn_pe = make_save_vcsn_year_mean('pe', np.nanmean, ['year', 'month'], outpath=os.path.join(outdir,'vcsn_pe_year_month.csv')).reset_index()
    vcsn_pe.loc[:, 'datetime'] = [pd.datetime(y, m, 15) for y, m, rch in vcsn_pe.itertuples(False, None)]
    vcsn_pe = vcsn_pe.set_index('datetime')
    peax.plot(vcsn_pe.index, vcsn_pe.pe, 'b', label='vcsn pe')

    rcp_past_pe = make_save_rcppast_year_amalg('pe', np.nanmean, ['year', 'month'],outpath=os.path.join(outdir,'rcp_past_pe_year_month.csv')).reset_index()
    rcp_past_pe.loc[:, 'datetime'] = [pd.datetime(y, m, 15) for y, m, rch in rcp_past_pe.itertuples(False, None)]
    rcp_past_pe = rcp_past_pe.set_index('datetime')
    peax.plot(rcp_past_pe.index, rcp_past_pe.pe, 'r', label='RCPpast pe')
    peax.legend()

    pefig, peax = plt.subplots()
    rcp_past_pe = make_save_rcppast_year_amalg('pe', np.nanmean, ['month'], outpath=os.path.join(outdir,'rcp_past_pe_month.csv'))
    vcsn_pe = make_save_vcsn_year_mean('pe', np.nanmean, ['month'],outpath=os.path.join(outdir,'vcsn_pe_month.csv'))

    peax.plot(rcp_past_pe.index, rcp_past_pe.pe, 'r', label='RCPpast pe')
    peax.plot(vcsn_pe.index, vcsn_pe.pe, 'b', label='vcsn pe')
    peax.legend()

    pefig, peax = plt.subplots()
    rcp_past_pe = make_save_rcppast_year_amalg('pe', np.nanmean, ['year'],outpath=os.path.join(outdir,'rcp_past_pe_year.csv'))
    vcsn_pe = make_save_vcsn_year_mean('pe', np.nanmean, ['year'], outpath=os.path.join(outdir,'vcsn_pe_year.csv'))

    peax.plot(rcp_past_pe.index, rcp_past_pe.pe, 'r', label='RCPpast pe')
    peax.plot(vcsn_pe.index, vcsn_pe.pe, 'b', label='vcsn pe')
    peax.legend()

    plt.show()
