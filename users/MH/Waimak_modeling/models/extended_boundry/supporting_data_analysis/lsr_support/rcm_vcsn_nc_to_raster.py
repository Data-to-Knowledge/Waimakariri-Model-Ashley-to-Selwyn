# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/10/2017 6:21 PM
"""

from __future__ import division
from core import env
import numpy as np
import netCDF4 as nc
import itertools
from glob import glob
import os
import datetime
import sys
import geopandas as gpd
from scipy import stats
import matplotlib.pyplot as plt
from copy import deepcopy
import statsmodels.formula.api as smf
import statsmodels.api as sm
import pandas as pd
import pickle

# dimensions (time, lat, lon) in a square to match the LRSM data zone make the dimensions of the grid identical
geom = gpd.read_file(
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\bound.shp"
).loc[0, 'geometry']
lat_bounds = geom.bounds[1], geom.bounds[3]
lon_bounds = geom.bounds[0], geom.bounds[2]
time_bounds = (datetime.datetime(1972, 1, 2), datetime.datetime(2005, 12, 31))
outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data"
if not os.path.exists(outdir):
    os.makedirs(outdir)


def make_month_year_rcm_data():
    paths = glob(env.gw_met_data(
        r"niwa_netcdf\climate_projections\RCPpast\*\PE_VCSN_*_RCPpast_1971_2005_south-island_p05_daily_ECan.nc"))
    for pe_path in paths:
        print(pe_path)
        model = os.path.basename(os.path.dirname(pe_path))

        # pull out data
        pe_data = nc.Dataset(pe_path)
        pe = np.array(pe_data.variables['pe'])
        pe[np.isclose(pe, pe_data.variables['pe']._FillValue)] = np.nan
        time = nc.num2date(np.array(pe_data.variables['time']), pe_data.variables['time'].units)
        lats, lons = np.array(pe_data['latitude']), np.array(pe_data['longitude'])

        # set bounds for lat lon and time
        lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
        lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))
        time_idx = np.where((time >= time_bounds[0]) & (time <= time_bounds[1]))

        # bound the data
        lats = lats[lat_idx]
        lons = lons[lon_idx]
        time = time[time_idx]

        pe = pe[time_idx[0], :, :][:, lat_idx[0], :][:, :, lon_idx[0]]

        years = [e.year for e in time]
        months = [e.month for e in time]
        iterator = list(itertools.product(set(years), set(months)))
        out_pe = np.zeros((len(iterator), len(lats), len(lons))) * np.nan

        for i, (year, month) in enumerate(iterator):
            # amalgimate the rcm to the month and year level
            amalg_idx = np.where(np.isclose(years, year) & np.isclose(months, month))
            temp = np.nanmean(pe[amalg_idx[0], :, :], 0)
            out_pe[i] = temp

        # save to netcdf4
        out_lats = lats
        out_lons = lons
        out_time = np.array([datetime.datetime(y, m, 15) for y, m in iterator])

        outpath = os.path.join(outdir, '{}_monthly_pe.nc'.format(model))
        _write_netcdf_year(outpath, out_lats, out_lons, out_time, out_pe, pe_path, model)


def _write_netcdf_year(outpath, out_lats, out_lons, out_time, out_pe, pe_path, model):
    outfile = nc.Dataset(outpath, 'w')

    # create dimensions
    outfile.createDimension('latitude', len(out_lats))
    outfile.createDimension('longitude', len(out_lons))
    outfile.createDimension('time', len(out_time))

    # create variables
    time = outfile.createVariable('time', 'f8', ('time',), fill_value=np.nan)
    time.setncatts({'units': 'days since 1972-01-02 00:00:00',
                    'long_name': 'time',
                    'missing_value': np.nan})
    time[:] = nc.date2num(out_time, time.units)

    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'latitude',
                   'missing_value': np.nan})
    lat[:] = out_lats

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'longitude',
                   'missing_value': np.nan})
    lon[:] = out_lons

    pe = outfile.createVariable('pe', 'f8', ('time', 'latitude', 'longitude'), fill_value=np.nan)
    pe.setncatts({'units': 'mmday^(-1)',
                  'long_name': 'Potential Evapotranspiration',
                  'missing_value': np.nan})
    pe[:] = out_pe

    # set global attributes
    outfile.description = 'calander month averaged pe for model: {}, for a limited area lon:{} , lat:{} and time: {}'.format(
        model, lon_bounds, lat_bounds, time_bounds)
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(pe_path, sys.argv[0])


def _write_netcdf_regress_seasonal(outpath, out_lats, out_lons, out_data, pe_path, model):
    outfile = nc.Dataset(outpath, 'w')

    # create dimensions
    outfile.createDimension('latitude', len(out_lats))
    outfile.createDimension('longitude', len(out_lons))
    outfile.createDimension('season', 4)

    # create variables
    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'latitude',
                   'missing_value': np.nan})
    lat[:] = out_lats

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'longitude',
                   'missing_value': np.nan})
    lon[:] = out_lons

    month = outfile.createVariable('season', 'i4', ('season'), fill_value=0)
    month.setncatts({'units': '0: djf, 1: mam, 2:jja, 3:son',
                     'long_name': 'month',
                     'missing_value': np.nan})
    month[:] = np.arange(4)

    for key, data in out_data.items():
        temp = outfile.createVariable(key, 'f8', ('season', 'latitude', 'longitude'), fill_value=np.nan)
        temp.setncatts({'units': 'unitless',
                        'long_name': key,
                        'missing_value': np.nan})
        temp[:] = data

    # set global attributes
    outfile.description = 'regression coeffients for vcsn and model : {}, for a limited area lon:{} , lat:{}'.format(
        model, lon_bounds, lat_bounds)
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(pe_path, sys.argv[0])


def _write_netcdf_regress_all(outpath, out_lats, out_lons, out_data, pe_path, model):
    outfile = nc.Dataset(outpath, 'w')

    # create dimensions
    outfile.createDimension('latitude', len(out_lats))
    outfile.createDimension('longitude', len(out_lons))

    # create variables
    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'latitude',
                   'missing_value': np.nan})
    lat[:] = out_lats

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lat.setncatts({'units': 'WGS84',
                   'long_name': 'longitude',
                   'missing_value': np.nan})
    lon[:] = out_lons

    for key, data in out_data.items():
        temp = outfile.createVariable(key, 'f8', ('latitude', 'longitude'), fill_value=np.nan)
        temp.setncatts({'units': 'unitless',
                        'long_name': key,
                        'missing_value': np.nan})
        temp[:] = data

    # set global attributes
    outfile.description = 'regression coeffients for vcsn and model : {}, for a limited area lon:{} , lat:{}'.format(
        model, lon_bounds, lat_bounds)
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(pe_path, sys.argv[0])


def make_month_year_vcsn_data():
    pe_path = r"Y:\VirtualClimate\vcsn_precip_et_2016-06-06.nc"
    print(pe_path)
    model = os.path.basename(os.path.dirname(pe_path))

    # pull out data
    pe_data = nc.Dataset(pe_path)
    pe = np.array(pe_data.variables['pe'])
    time = nc.num2date(np.array(pe_data.variables['time']), pe_data.variables['time'].units)
    lats, lons = np.array(pe_data['latitude']), np.array(pe_data['longitude'])
    lats = np.flip(lats, 0)
    pe = np.swapaxes(pe, 0, 2)
    pe = np.flip(pe, axis=1)

    # set bounds for lat lon and time
    lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
    lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))
    time_idx = np.where((time >= time_bounds[0]) & (time <= time_bounds[1]))

    # bound the data
    lats = lats[lat_idx]
    lons = lons[lon_idx]
    time = time[time_idx]

    pe = pe[time_idx[0], :, :][:, lat_idx[0], :][:, :, lon_idx[0]]

    years = [e.year for e in time]
    months = [e.month for e in time]
    iterator = list(itertools.product(set(years), set(months)))
    out_pe = np.zeros((len(iterator), len(lats), len(lons))) * np.nan

    for i, (year, month) in enumerate(iterator):
        # amalgimate the rcm to the month and year level
        amalg_idx = np.where(np.isclose(years, year) & np.isclose(months, month))
        temp = np.nanmean(pe[amalg_idx[0], :, :], 0)
        out_pe[i] = temp
    # apply course step change correction for the added inland site in 2000  see jen's memo: Comparison of rainfall and
    # evaporation data with the Virtual Climate Station network (VCSN) in the Ashley-Waimakariri zone
    out_time = np.array([datetime.datetime(y, m, 15) for y, m in iterator])
    idx_2000 = out_time < datetime.datetime(2001, 1, 1)
    for i, j in itertools.product(range(out_pe.shape[1]), range(out_pe.shape[2])):
        temp = pd.DataFrame(index=out_time, data=out_pe[:, i, j], columns=['pe'])
        temp = temp.resample('A').sum()
        scale = temp.loc[temp.index >= datetime.datetime(2001, 1, 1), 'pe'].mean() / temp.loc[
            temp.index < datetime.datetime(2001, 1, 1), 'pe'].mean()
        out_pe[idx_2000, i, j] *= scale

    # save to netcdf4
    out_lats = lats
    out_lons = lons

    outpath = os.path.join(outdir, '{}_monthly_pe.nc'.format(model))
    _write_netcdf_year(outpath, out_lats, out_lons, out_time, out_pe, pe_path, model)


def make_slope_intercept_grid_seasonal():
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\*-*.nc")
    vcsn_path = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\VirtualClimate_monthly_pe.nc"
    vcsn_data = nc.Dataset(vcsn_path)
    vcsn_pe = np.array(vcsn_data.variables['pe'])
    lats = np.array(vcsn_data.variables['latitude'])
    lons = np.array(vcsn_data.variables['longitude'])
    months = np.array([e.month for e in nc.num2date(np.array(vcsn_data['time']), vcsn_data['time'].units)])
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\seasonal"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for path in paths:
        print(path)
        model = os.path.basename(path).replace('_monthly_pe.nc', '')
        outpath = os.path.join(outdir, '{}_vcsn_regression.nc'.format(model))
        rcm_data = nc.Dataset(path)
        rcm_pe = np.array(rcm_data.variables['pe'])

        if rcm_pe.shape != vcsn_pe.shape:
            raise ValueError('rcm and vcsn shapes do not match')

        outshape = (4, rcm_pe.shape[1], rcm_pe.shape[2])
        counts = np.zeros(outshape) * np.nan
        slopes = np.zeros(outshape) * np.nan
        intercepts = np.zeros(outshape) * np.nan
        r_values = np.zeros(outshape) * np.nan
        p_values = np.zeros(outshape) * np.nan
        std_errs = np.zeros(outshape) * np.nan
        seasons = {0: (12, 1, 2), 1: (3, 4, 5), 2: (6, 7, 8), 3: (9, 10, 11)}
        for month, latidx, lonidx in itertools.product(range(4), range(rcm_pe.shape[1]), range(rcm_pe.shape[2])):
            season_months = seasons[month]
            x = rcm_pe[:, latidx, lonidx]
            y = vcsn_pe[:, latidx, lonidx]
            mask = np.isfinite(x) & np.isfinite(y) & np.in1d(months, season_months)
            count = mask.sum()
            if count == 0:
                slope, intercept, r_value, p_value, std_err = np.nan, np.nan, np.nan, np.nan, np.nan
            else:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask], y[mask])
            counts[month, latidx, lonidx] = count
            slopes[month, latidx, lonidx] = slope
            intercepts[month, latidx, lonidx] = intercept
            r_values[month, latidx, lonidx] = r_value
            p_values[month, latidx, lonidx] = p_value
            std_errs[month, latidx, lonidx] = std_err

        outdata = {'slope': slopes, 'intercept': intercepts, 'r_value': r_values,
                   'p_value': p_values, 'std_err': std_errs, 'count': counts}
        _write_netcdf_regress_seasonal(outpath, lats, lons, outdata, path, model)


def make_slope_intercept_grid_all():
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\*-*.nc")
    vcsn_path = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\VirtualClimate_monthly_pe.nc"
    vcsn_data = nc.Dataset(vcsn_path)
    vcsn_pe = np.array(vcsn_data.variables['pe'])
    lats = np.array(vcsn_data.variables['latitude'])
    lons = np.array(vcsn_data.variables['longitude'])
    months = np.array([e.month for e in nc.num2date(np.array(vcsn_data['time']), vcsn_data['time'].units)])
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\all"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for path in paths:
        print(path)
        model = os.path.basename(path).replace('_monthly_pe.nc', '')
        outpath = os.path.join(outdir, '{}_vcsn_regression.nc'.format(model))
        rcm_data = nc.Dataset(path)
        rcm_pe = np.array(rcm_data.variables['pe'])

        if rcm_pe.shape != vcsn_pe.shape:
            raise ValueError('rcm and vcsn shapes do not match')

        outshape = (rcm_pe.shape[1], rcm_pe.shape[2])
        counts = np.zeros(outshape) * np.nan
        slopes = np.zeros(outshape) * np.nan
        intercepts = np.zeros(outshape) * np.nan
        r_values = np.zeros(outshape) * np.nan
        p_values = np.zeros(outshape) * np.nan
        std_errs = np.zeros(outshape) * np.nan
        rcm_wp = np.zeros(outshape) * np.nan
        vcsn_wp = np.zeros(outshape) * np.nan
        for latidx, lonidx in itertools.product(range(rcm_pe.shape[1]), range(rcm_pe.shape[2])):
            x = np.log(rcm_pe[:, latidx, lonidx])
            y = np.log(vcsn_pe[:, latidx, lonidx])
            mask = np.isfinite(x) & np.isfinite(y)
            count = mask.sum()
            if count == 0:
                slope, intercept, r_value, p_value, std_err = np.nan, np.nan, np.nan, np.nan, np.nan
            else:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask], y[mask])
                rcm_w, rcm_wpv = stats.shapiro(x[mask])
                vcsn_w, vcsn_wpv = stats.shapiro(y[mask])

            counts[latidx, lonidx] = count
            slopes[latidx, lonidx] = slope
            intercepts[latidx, lonidx] = intercept
            r_values[latidx, lonidx] = r_value
            p_values[latidx, lonidx] = p_value
            std_errs[latidx, lonidx] = std_err
            rcm_wp[latidx, lonidx] = rcm_wpv
            vcsn_wp[latidx, lonidx] = vcsn_wpv

        outdata = {'slope': slopes, 'intercept': intercepts, 'r_value': r_values,
                   'p_value': p_values, 'std_err': std_errs, 'count': counts, 'shapiro_rcm': rcm_wp,
                   'shapiro_vcsn': vcsn_wp}
        _write_netcdf_regress_all(outpath, lats, lons, outdata, path, model)


def make_slope_intercept_grid_all_ols():
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\*-*.nc")
    vcsn_path = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\VirtualClimate_monthly_pe.nc"
    vcsn_data = nc.Dataset(vcsn_path)
    vcsn_pe = np.array(vcsn_data.variables['pe'])
    lats = np.array(vcsn_data.variables['latitude'])
    lons = np.array(vcsn_data.variables['longitude'])
    months = np.array([e.month for e in nc.num2date(np.array(vcsn_data['time']), vcsn_data['time'].units)])
    outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\all_ols"
    month_to_season = {12: 0, 1: 0, 2: 0,  # summer
                       3: 1, 4: 1, 5: 1,  # spring
                       6: 2, 7: 2, 8: 2,  # fall
                       9: 3, 10: 3, 11: 3}  # winter
    season_names = {0: 'summer', 1: 'fall', 2: 'winter', 3: 'spring'}
    season = np.array([month_to_season[e] for e in months])
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for path in paths:
        print(path)
        model = os.path.basename(path).replace('_monthly_pe.nc', '')
        outpath = os.path.join(outdir, '{}_vcsn_regression.nc'.format(model))
        rcm_data = nc.Dataset(path)
        rcm_pe = np.array(rcm_data.variables['pe'])

        if rcm_pe.shape != vcsn_pe.shape:
            raise ValueError('rcm and vcsn shapes do not match')

        outshape = (rcm_pe.shape[1], rcm_pe.shape[2])
        counts = np.zeros(outshape) * np.nan
        slopes = np.zeros(outshape) * np.nan
        intercepts = np.zeros(outshape) * np.nan
        t1s = np.zeros(outshape) * np.nan
        t2s = np.zeros(outshape) * np.nan
        t3s = np.zeros(outshape) * np.nan
        r_values = np.zeros(outshape) * np.nan
        fp_values = np.zeros(outshape) * np.nan
        p_values_x = np.zeros(outshape) * np.nan
        p_values_t1 = np.zeros(outshape) * np.nan
        p_values_t2 = np.zeros(outshape) * np.nan
        p_values_t3 = np.zeros(outshape) * np.nan
        p_values_intercept = np.zeros(outshape) * np.nan
        regressions = np.zeros(outshape, object)
        for latidx, lonidx in itertools.product(range(rcm_pe.shape[1]), range(rcm_pe.shape[2])):
            x = np.log(rcm_pe[:, latidx, lonidx])
            y = np.log(vcsn_pe[:, latidx, lonidx])
            mask = np.isfinite(x) & np.isfinite(y)
            count = mask.sum()
            if count == 0:
                slope, intercept, r_value, = np.nan, np.nan, np.nan
                t1, t2, t3, p_x, p_t1, p_t2, p_t3, p_intercept, fpv = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                res = None
            else:
                temp_data = pd.DataFrame({'x': x, 'y': y, 'season': season})
                mod = smf.ols(formula='y ~ x + C(season)', data=temp_data)
                res = mod.fit()
                slope = res.params.loc['x']
                intercept = res.params.loc['Intercept']
                t1 = res.params.loc['C(season)[T.1]']
                t2 = res.params.loc['C(season)[T.2]']
                t3 = res.params.loc['C(season)[T.3]']
                r_value = res.rsquared_adj

                p_x = res.pvalues.loc['x']
                p_t1 = res.pvalues.loc['C(season)[T.1]']
                p_t2 = res.pvalues.loc['C(season)[T.2]']
                p_t3 = res.pvalues.loc['C(season)[T.3]']
                p_intercept = res.pvalues.loc['Intercept']
                fpv = res.f_pvalue
            counts[latidx, lonidx] = count
            slopes[latidx, lonidx] = slope
            t1s[latidx, lonidx] = t1
            t2s[latidx, lonidx] = t2
            t3s[latidx, lonidx] = t3
            intercepts[latidx, lonidx] = intercept
            r_values[latidx, lonidx] = r_value
            p_values_x[latidx, lonidx] = p_x
            p_values_t1[latidx, lonidx] = p_t1
            p_values_t2[latidx, lonidx] = p_t2
            p_values_t3[latidx, lonidx] = p_t3
            p_values_intercept[latidx, lonidx] = p_intercept
            fp_values[latidx, lonidx] = fpv
            regressions[latidx, lonidx] = res

        variables = ['counts', 'slopes', 't1s', 't2s', 't3s', 'intercepts', 'r_values', 'p_values_x', 'p_values_t1',
                     'p_values_t2', 'p_values_t3', 'p_values_intercept', 'fp_values']

        outdata = {}
        for e in variables:
            outdata[e] = eval(e)

        _write_netcdf_regress_all(outpath, lats, lons, outdata, path, model)
        pickle.dump(regressions, open(os.path.join(outdir, 'regressions_{}.p'.format(model)), 'w'))


def visualise_regress_data_seasonal(path, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    dataset = nc.Dataset(path)
    print(os.path.basename(path))
    print(set(np.array(dataset.variables['count']).flatten()))
    seasons = {0: 'djf', 1: 'mam', 2: 'jja', 3: 'son'}
    for var in ['slope', 'intercept', 'r_value']:
        data = np.array(dataset.variables[var])
        fig, axs = plt.subplots(2, 2, figsize=(18.5, 9.5))
        for i, ax in enumerate(axs.flatten()):
            temp = ax.pcolormesh(np.flipud(data[i]), cmap='RdBu', vmin=np.nanmin(data), vmax=np.nanmax(data))
            ax.set_aspect('equal')
            ax.set_title(seasons[i])
        fig.colorbar(temp)
        model = os.path.basename(path).split('_')[0]
        fig.suptitle('{} {}'.format(model, var))
        fig.savefig(os.path.join(outdir, '{}_{}.png'.format(model, var)))


def visualise_regress_data_all(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\all\*.nc")

    for var in ['shapiro_rcm', 'shapiro_vcsn', 'slope', 'intercept', 'r_value', 'std_err']:
        maxes, mins = [], []
        for path in paths:
            dataset = nc.Dataset(path)
            data = np.array(dataset.variables[var])
            maxes.append(np.nanmax(data))
            mins.append(np.nanmin(data))

        fig, axs = plt.subplots(2, 3, figsize=(18.5, 9.5))
        for i, ax, path in zip(range(len(paths)), axs.flatten(), paths):
            dataset = nc.Dataset(path)
            data = np.array(dataset.variables[var])
            temp = ax.pcolormesh(np.flipud(data), cmap='RdBu', vmin=np.nanmin(mins), vmax=np.nanmax(maxes))
            ax.set_aspect('equal')
            ax.set_title(os.path.basename(path).split('_')[0])
        fig.colorbar(temp)
        fig.suptitle(var)
        fig.savefig(os.path.join(outdir, 'all_{}.png'.format(var)))


def visualise_regress_data_all_ols(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\all_ols\*.nc")

    for var in ['counts', 'slopes', 't1s', 't2s', 't3s', 'intercepts', 'r_values', 'p_values_x', 'p_values_t1',
                'p_values_t2', 'p_values_t3', 'p_values_intercept', 'fp_values']:
        maxes, mins = [], []
        for path in paths:
            dataset = nc.Dataset(path)
            data = np.array(dataset.variables[var])
            maxes.append(np.nanmax(data))
            mins.append(np.nanmin(data))

        fig, axs = plt.subplots(2, 3, figsize=(18.5, 9.5))
        for i, ax, path in zip(range(len(paths)), axs.flatten(), paths):
            dataset = nc.Dataset(path)
            data = np.array(dataset.variables[var])
            temp = ax.pcolormesh(np.flipud(data), cmap='RdBu', vmin=np.nanmin(mins), vmax=np.nanmax(maxes))
            ax.set_aspect('equal')
            ax.set_title(os.path.basename(path).split('_')[0])
        fig.colorbar(temp)
        fig.suptitle(var)
        fig.savefig(os.path.join(outdir, 'all_{}.png'.format(var)))


def plot_example_sites(outdir, log_tranform=False):  # log is better
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    shp_sites = gpd.read_file(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\sites.shp")
    sites = {name: (geom.x, geom.y) for name, geom in
             shp_sites.loc[:, ['site_name', 'geometry']].itertuples(False, None)}
    data = {}
    paths = glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\*-*.nc")
    for path in paths:
        temp = nc.Dataset(path)
        model = os.path.basename(path).replace('_monthly_pe.nc', '')
        data[model] = np.array(temp.variables['pe'])

    base_cols = ['c', 'y', 'm', 'k', 'g', 'b']
    colors = {key: col for key, col in zip(data.keys(), base_cols)}

    vcsn_path = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\VirtualClimate_monthly_pe.nc"
    vcsn_data = nc.Dataset(vcsn_path)
    vcsn = np.array(vcsn_data.variables['pe'])
    lats = np.array(vcsn_data.variables['latitude'])
    lons = np.array(vcsn_data.variables['longitude'])
    time = nc.num2date(np.array(vcsn_data.variables['time']), vcsn_data.variables['time'].units)
    month_to_season = {12: 0, 1: 0, 2: 0,  # summer
                       3: 1, 4: 1, 5: 1,  # spring
                       6: 2, 7: 2, 8: 2,  # fall
                       9: 3, 10: 3, 11: 3}  # winter
    season_colors = {0: 'y', 1: 'r', 2: 'b', 3: 'g'}
    season_names = {0: 'summer', 1: 'fall', 2: 'winter', 3: 'spring'}
    season = np.array([month_to_season[e.month] for e in time])

    for site, (lon, lat) in sites.items():
        lat_idx = np.argmin(np.abs(np.abs(lats) - np.abs(lat)))
        lon_idx = np.argmin(np.abs(np.abs(lons) - np.abs(lon)))
        print(lat, lat_idx, lon, lon_idx)
        vcs_for_site = vcsn[:, lat_idx, lon_idx]
        fig, axts = plt.subplots(figsize=(18.5, 9.5))
        fig2, axs = plt.subplots(2, 3, figsize=(18.5, 9.5))
        for ax, (m, dat) in zip(axs.flatten(), data.items()):
            rcm_site = deepcopy(dat[:, lat_idx, lon_idx])
            # make a ts plot
            axts.plot(time, rcm_site, c=colors[m], label=m)
            axts.plot(time, vcs_for_site, c='r', label='vcsn')
            axts.legend()
            axts.set_ylabel('pet mm/day')
            axts.set_xlabel('time')

            # make a linear regression plot
            if log_tranform:
                rcm_site = np.log(rcm_site)
                vcsn_for_site_use = np.log(vcs_for_site)
                des = 'log'
            else:
                vcsn_for_site_use = vcs_for_site
                des = 'non-log'
            slope, intercept, r_value, p_value, std_err = stats.linregress(rcm_site, vcsn_for_site_use)
            for i in range(4):
                idx = season == i
                ax.scatter(rcm_site[idx], vcsn_for_site_use[idx], c=season_colors[i], marker='.', label=season_names[i])
            abline(slope, intercept, ax)
            ax.set_title(m)
            ax.set_ylabel = 'vcsn'
            ax.set_xlabel = 'rcm'
            ax.annotate('y={:.2f}x+{:.2f}'.format(slope, intercept), xy=(0.05, 0.95), xycoords='axes fraction')
        fig.suptitle(site)
        fig.savefig(os.path.join(outdir, 'ts_{}__plot.png'.format(site)))
        handles, labels = fig2.axes[0].get_legend_handles_labels()
        fig2.legend(handles, labels, 'upper right')
        fig2.suptitle('{}-{}'.format(des, site))
        fig2.savefig(os.path.join(outdir, 'reg_{}_plot.png'.format(site)))


def abline(slope, intercept, ax, color=None):
    """Plot a line from slope and intercept"""
    x_vals = np.array(ax.get_xlim())
    y_vals = intercept + slope * x_vals
    ax.plot(x_vals, y_vals, '--', c=color)


if __name__ == '__main__':
    # make_month_year_rcm_data()
    make_month_year_vcsn_data()
    make_slope_intercept_grid_all_ols()
    plot_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\plots"
    # plot all
    visualise_regress_data_all_ols(os.path.join(plot_dir, 'all_ols'))

    plt_regress = True
    if plt_regress:
        plot_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\plots"
        # plot all
        visualise_regress_data_all(os.path.join(plot_dir, 'all'))

        # plot seasonal
        paths = glob(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\seasonal\*.nc")
        for path in paths:
            visualise_regress_data_seasonal(path, os.path.join(plot_dir, 'seasonal'))
            # plot examples sites
        plot_example_sites(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\plots\notlog_sites")
        plot_example_sites(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\plots\log_sites",
            True)
