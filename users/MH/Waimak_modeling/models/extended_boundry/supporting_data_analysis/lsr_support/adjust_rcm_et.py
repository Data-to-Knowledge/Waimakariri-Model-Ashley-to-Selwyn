# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/10/2017 3:09 PM
"""

from __future__ import division
from core import env
import geopandas as gpd
import numpy as np
from glob import glob
import os
import netCDF4 as nc
import itertools
import pickle
import statsmodels
import datetime
import sys
import matplotlib.pyplot as plt
import pandas as pd

geom = gpd.read_file(
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\bound.shp"
).loc[0, 'geometry']
lat_bounds = geom.bounds[1], geom.bounds[3]
lon_bounds = geom.bounds[0], geom.bounds[2]


def adjust_rcm_pet(baseoutdir):
    """
    adjust the rcm_pe by regression calculated in rcm_cvsn_nc_to_raster
    :param outdir:
    :return:
    """
    # get the path for the rcm data and load
    paths = glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCP8.5\*\PE*.nc"))
    paths.extend(glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCP4.5\*\PE*.nc")))
    paths.extend(glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCPpast\*\PE*.nc")))

    for pe_path in paths:
        model = os.path.basename(os.path.dirname(pe_path))
        rcm = os.path.basename(os.path.dirname((os.path.dirname(pe_path))))
        outdir = os.path.join(baseoutdir,rcm,model)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        print(pe_path)

        # pull out data
        pe_data = nc.Dataset(pe_path)
        pe = np.array(pe_data.variables['pe'])
        pe[np.isclose(pe, pe_data.variables['pe']._FillValue)] = np.nan
        lats, lons = np.array(pe_data['latitude']), np.array(pe_data['longitude'])
        time = nc.num2date(np.array(pe_data.variables['time']), pe_data.variables['time'].units)
        month_to_season = {12: 0, 1: 0, 2: 0,  # summer
                           3: 1, 4: 1, 5: 1,  # spring
                           6: 2, 7: 2, 8: 2,  # fall
                           9: 3, 10: 3, 11: 3}  # winter
        season_names = {0: 'summer', 1: 'fall', 2: 'winter', 3: 'spring'}
        seasons = np.array([month_to_season[e.month] for e in time])

        # set bounds for lat lon
        lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
        lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))

        # bound the data
        lats = lats[lat_idx]
        lons = lons[lon_idx]
        pe = pe[:, lat_idx[0], :][:, :, lon_idx[0]]

        # get regression data
        reg_path = env.sci(
            "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model simulations and results/ex_bd_va/niwa_data_explore/montly_nc_data/regressions/all_ols/regressions_{}.p".format(
                model))
        reg_data = pickle.load(open(reg_path))
        # and apply correction to the data
        new_pe_data = np.zeros(pe.shape) * np.nan
        for latidx, lonidx in itertools.product(range(pe.shape[1]), range(pe.shape[2])):
            reg = reg_data[latidx, lonidx]
            if reg is None:
                continue
            new_pe_data[:, latidx, lonidx] = predict(reg, pe[:, latidx, lonidx], seasons)
        outpath = os.path.join(outdir, 'PE_waimak_adjusted_{}.nc'.format(model))

        # save as a new netcdf similar to the old netcdf
        _write_netcdf_year(outpath, lats, lons, time, new_pe_data, pe_path=pe_path, model=model)



def _write_netcdf_year(outpath, out_lats, out_lons, out_time, out_pe, pe_path, model):
    outfile = nc.Dataset(outpath, 'w')

    # create dimensions
    outfile.createDimension('latitude', len(out_lats))
    outfile.createDimension('longitude', len(out_lons))
    outfile.createDimension('time', len(out_time))

    # create variables
    time = outfile.createVariable('time', 'f8', ('time',), fill_value=np.nan)
    time.setncatts({'units': 'days since 1970-01-01 00:00:00',
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
    outfile.description = ('corrected pe for model: {}, for a limited area lon:{} , lat:{}\n'
                          ' Correction based on a model and grid cell specific OLS regression between the RCM (x) and '
                          'VCSN (y) monthly averaged pe data from 1972 to 2006 and the catagorised season '
                           '(summer = djf, fall = mam, winter = jja, spring: son)'.format(model, lon_bounds,
                                                                                          lat_bounds))
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(pe_path, sys.argv[0])


def _write_netcdf_rain(outpath, out_lats, out_lons, out_time, out_rain, rain_path, model):
    outfile = nc.Dataset(outpath, 'w')

    # create dimensions
    outfile.createDimension('latitude', len(out_lats))
    outfile.createDimension('longitude', len(out_lons))
    outfile.createDimension('time', len(out_time))

    # create variables
    time = outfile.createVariable('time', 'f8', ('time',), fill_value=np.nan)
    time.setncatts({'units': 'days since 1970-01-01 00:00:00',
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

    rain = outfile.createVariable('rain', 'f8', ('time', 'latitude', 'longitude'), fill_value=np.nan)
    rain.setncatts({'units': 'kg m-2',
                  'long_name': 'precipitation_amount (mm)',
                  'missing_value': np.nan,
                  'description': 'total value over previous 24 hours'})
    rain[:] = out_rain

    # set global attributes
    outfile.description = ('copied rainfall data for model: {}, for a limited area lon:{} , lat:{}'.format(model,
                                                                                                             lon_bounds,
                                                                                                            lat_bounds))
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'original data: {}, script: {}'.format(rain_path, sys.argv[0])


def predict(regression, pe_data, season):
    """
    if a wrapper to run the regression
    :param regression: statsmodels.regression.linear_model.RegressionResultsWrapper
    :param pe_data: (time,)
    :param season: (time, )
    :return:
    """
    if not isinstance(regression, statsmodels.regression.linear_model.RegressionResultsWrapper):
        raise ValueError('regression must be a regression results wrapper')
    # data was log/log transform
    if season.shape != pe_data.shape:
        raise ValueError('season shape must match shape the pe_data shape')
    pe_data = np.log(pe_data)
    temp = np.e ** regression.predict({'x': pe_data, 'season': season})

    return temp

def copy_area_precip_data(baseoutdir):
    # get the path for the rcm data and load
    paths = glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCP8.5\*\TotalPrecipCorr*.nc"))
    paths.extend(glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCP4.5\*\TotalPrecipCorr*.nc")))
    paths.extend(glob(env.gw_met_data(r"niwa_netcdf\climate_projections\RCPpast\*\TotalPrecipCorr*.nc")))
    for rain_path in paths:
        rcm = os.path.basename(os.path.dirname((os.path.dirname(rain_path))))
        model = os.path.basename(os.path.dirname(rain_path))
        outdir = os.path.join(baseoutdir,rcm,model)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        print(rain_path)

        # pull out data
        rain_data = nc.Dataset(rain_path)
        rain = np.array(rain_data.variables['rain'])
        rain[np.isclose(rain, rain_data.variables['rain']._FillValue)] = np.nan
        lats, lons = np.array(rain_data['latitude']), np.array(rain_data['longitude'])
        time = nc.num2date(np.array(rain_data.variables['time']), rain_data.variables['time'].units)

        # set bounds for lat lon
        lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
        lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))

        # bound the data
        lats = lats[lat_idx]
        lons = lons[lon_idx]
        rain = rain[:, lat_idx[0], :][:, :, lon_idx[0]]

        outpath = os.path.join(outdir, 'TotalPrecipCorr_waimak_adjusted_{}.nc'.format(model))

        # save as a new netcdf similar to the old netcdf
        _write_netcdf_rain(outpath, lats, lons, time, rain, rain_path=rain_path, model=model)



def check_new_rcm_data(base_dir):
    # monthly and annual averaged ET time periods at sites all rcms
    outdir = os.path.join(base_dir,'site_plots')
    rcps = ['RCPpast','RCP4.5','RCP8.5']
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    shp_sites = gpd.read_file(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\sites.shp")
    sites = {name: (geom.x, geom.y) for name, geom in
             shp_sites.loc[:, ['site_name', 'geometry']].itertuples(False, None)}
    for site, (lon, lat) in sites.items():
        print(site)
        for style in ['year','month','10year']:
            print(style)
            fig, axs = plt.subplots(nrows=3,figsize=(18.5, 9.5))
            for rcp,ax in zip(rcps, axs.flatten()):
                data = {}
                all_times = {}
                paths = glob(
                    r"K:\niwa_netcdf\rain_waimak_corrected_pe_data\{}\*\PE*.nc".format(rcp))
                for path in paths:
                    temp = nc.Dataset(path)
                    model = os.path.basename(path).replace('_monthly_pe.nc', '')
                    data[model] = np.array(temp.variables['pe'])
                    all_times[model] = nc.num2date(np.array(temp.variables['time']), temp.variables['time'].units)

                base_cols = ['c', 'y', 'm', 'k', 'g', 'b']
                colors = {key: col for key, col in zip(data.keys(), base_cols)}
                colors['vcsn'] = 'r'
                if rcp == 'RCPpast':
                    pe_path = r"Y:\VirtualClimate\vcsn_precip_et_2016-06-06.nc"
                    pe_data = nc.Dataset(pe_path)
                    pe = np.array(pe_data.variables['pe'])
                    vcsn_time = nc.num2date(np.array(pe_data.variables['time']), pe_data.variables['time'].units)
                    lats, lons = np.array(pe_data['latitude']), np.array(pe_data['longitude'])
                    lats = np.flip(lats, 0)
                    pe = np.swapaxes(pe, 0, 2)
                    pe = np.flip(pe, axis=1)

                    # set bounds for lat lon and time
                    lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
                    lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))

                    # bound the data
                    lats = lats[lat_idx]
                    lons = lons[lon_idx]

                    pe = pe[:, lat_idx[0], :][:, :, lon_idx[0]]
                    data['vcsn'] = pe
                    all_times['vcsn'] = vcsn_time

                lats = np.array(temp.variables['latitude'])
                lons = np.array(temp.variables['longitude'])

                lat_idx = np.argmin(np.abs(np.abs(lats) - np.abs(lat)))
                lon_idx = np.argmin(np.abs(np.abs(lons) - np.abs(lon)))
                print(lat, lat_idx, lon, lon_idx)
                for m, dat in data.items():
                    rcm_site_df = pd.DataFrame(index=all_times[m], data=dat[:, lat_idx, lon_idx], columns=['pet'])
                    if style == 'month':
                        plt_data = rcm_site_df.resample('M').mean()
                    elif style == 'year':
                        plt_data = rcm_site_df.resample('A').mean()
                    elif style == '10year':
                        plt_data = rcm_site_df.resample('10A').mean()
                    # make a ts plot
                    ax.plot(plt_data.index, plt_data.pet, c=colors[m], label=m)
                    ax.set_ylabel('pet mm/day')
                    ax.set_xlabel('time')
                    ax.set_title(rcp)

            fig.suptitle(site)
            handles, labels = fig.axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, 'upper right')
            fig.savefig(os.path.join(outdir, 'mean_PE_{}_ts_{}__plot.png'.format(style, site)))



def check_new_rcm_data_rain(base_dir):
    # monthly and annual averaged ET time periods at sites all rcms
    outdir = os.path.join(base_dir,'site_plots')
    rcps = ['RCPpast','RCP4.5','RCP8.5']
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    shp_sites = gpd.read_file(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\niwa_data_explore\montly_nc_data\regressions\sites.shp")
    sites = {name: (geom.x, geom.y) for name, geom in
             shp_sites.loc[:, ['site_name', 'geometry']].itertuples(False, None)}
    for site, (lon, lat) in sites.items():
        print(site)
        for style in ['year','month','10year']:
            print(style)
            fig, axs = plt.subplots(nrows=3,figsize=(18.5, 9.5))
            for rcp,ax in zip(rcps, axs.flatten()):
                data = {}
                all_times = {}
                paths = glob(
                    r"K:\niwa_netcdf\rain_waimak_corrected_pe_data\{}\*\Tot*.nc".format(rcp))
                for path in paths:
                    temp = nc.Dataset(path)
                    model = os.path.basename(os.path.dirname(path))
                    data[model] = np.array(temp.variables['rain'])
                    all_times[model] = nc.num2date(np.array(temp.variables['time']), temp.variables['time'].units)

                base_cols = ['c', 'y', 'm', 'k', 'g', 'b']
                colors = {key: col for key, col in zip(data.keys(), base_cols)}
                colors['vcsn'] = 'r'
                if rcp == 'RCPpast':
                    pe_path = r"Y:\VirtualClimate\vcsn_precip_et_2016-06-06.nc"
                    pe_data = nc.Dataset(pe_path)
                    pe = np.array(pe_data.variables['rain'])
                    vcsn_time = nc.num2date(np.array(pe_data.variables['time']), pe_data.variables['time'].units)
                    lats, lons = np.array(pe_data['latitude']), np.array(pe_data['longitude'])
                    lats = np.flip(lats, 0)
                    pe = np.swapaxes(pe, 0, 2)
                    pe = np.flip(pe, axis=1)

                    # set bounds for lat lon and time
                    lat_idx = np.where((lats >= lat_bounds[0]) & (lats <= lat_bounds[1]))
                    lon_idx = np.where((lons >= lon_bounds[0]) & (lons <= lon_bounds[1]))

                    # bound the data
                    lats = lats[lat_idx]
                    lons = lons[lon_idx]

                    pe = pe[:, lat_idx[0], :][:, :, lon_idx[0]]
                    data['vcsn'] = pe
                    all_times['vcsn'] = vcsn_time

                lats = np.array(temp.variables['latitude'])
                lons = np.array(temp.variables['longitude'])

                lat_idx = np.argmin(np.abs(np.abs(lats) - np.abs(lat)))
                lon_idx = np.argmin(np.abs(np.abs(lons) - np.abs(lon)))
                print(lat, lat_idx, lon, lon_idx)
                for m, dat in data.items():
                    rcm_site_df = pd.DataFrame(index=all_times[m], data=dat[:, lat_idx, lon_idx], columns=['pet'])
                    if style == 'month':
                        plt_data = rcm_site_df.resample('M').sum()
                    elif style == 'year':
                        plt_data = rcm_site_df.resample('A').sum()
                    elif style == '10year':
                        plt_data = rcm_site_df.resample('10A').mean()
                    # make a ts plot
                    ax.plot(plt_data.index, plt_data.pet, c=colors[m], label=m)
                    ax.set_ylabel('rain mm/day')
                    ax.set_xlabel('time')
                    ax.set_title(rcp)

            fig.suptitle(site)
            handles, labels = fig.axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, 'upper right')
            fig.savefig(os.path.join(outdir, 'rain_sum_{}_ts_{}__plot.png'.format(style, site)))





if __name__ == '__main__':
    #copy_area_precip_data(env.gw_met_data("niwa_netcdf/rain_waimak_corrected_pe_data"))
    adjust_rcm_pet(env.gw_met_data("niwa_netcdf/rain_waimak_corrected_pe_data"))
    check_new_rcm_data(r"K:\niwa_netcdf\rain_waimak_corrected_pe_data\check_plots")
    #check_new_rcm_data_rain(r"K:\niwa_netcdf\rain_waimak_corrected_pe_data\check_plots")
