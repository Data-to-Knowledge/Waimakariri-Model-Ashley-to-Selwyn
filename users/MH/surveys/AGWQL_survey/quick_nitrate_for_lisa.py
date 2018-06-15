# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 24/05/2018 9:38 AM
"""

from __future__ import division
from core import env
import pandas as pd
import matplotlib.pyplot as plt
from core.stats.mann_kendall import mann_kendall_obj
from plotting import plot_background
import os
import matplotlib.patches as patches
import numpy as np


def _get_first(x):
    return x.iloc[0]


def make_n_trends(data, outdir, site_col, datetime_col, value_col, easting_col, northing_col,
                  gis_dir="P:/Groundwater/Annual groundwater quality survey 2016_true_copy/GIS_data"):
    """

    :param data: pandas data frame
    :param outdir: directory to send evertying
    :param site_col: the column with the site identifiers
    :param datetime_col: the column with the datetime object
    :param valuecol: the column with the values
    :param easting_col: the col with the lon
    :param northing_col: the col with the lat
    :param gis_dir: directory with the shapefiles in it
    :return:
    """
    sites = data[site_col].unique()
    # plot data
    plot_dir = os.path.join(outdir, 'plots')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    for site in sites:
        temp = data.loc[data[site_col] == site].set_index(datetime_col)
        fig, ax = plt.subplots()
        ax.scatter(temp.index, temp[value_col])
        ax.set_title(site)
        ax.set_ylabel('NO3-N')
        ax.set_xlabel('time')
        fig.savefig(os.path.join(plot_dir, site.replace('/', '_') + '.png'))
        plt.close(fig)

    # make trends
    trend_data = pd.DataFrame(index=sites, columns=['trend', 'count', 'northing', 'easting'])
    trend_data.loc[:, 'count'] = data.groupby(site_col).count()
    temp = data.groupby(site_col).aggregate({easting_col: _get_first, northing_col: _get_first})
    trend_data.loc[:, 'lat'] = temp[northing_col]
    trend_data.loc[:, 'lon'] = temp[easting_col]
    for site in sites:
        temp = mann_kendall_obj(data.loc[data[site_col]==site].set_index(datetime_col), data_col=value_col)
        trend_data.loc[site, 'trend'] = temp.trend
    trend_data.to_csv(os.path.join(outdir,'n_trends.csv'))

    # plot overview
    fig, ax = plot_background(gis_dir)
    tempdata = trend_data[trend_data['trend'] == 'no trend']
    plt.plot(tempdata.lon, tempdata.lat, color='grey', marker='o', label='No Trend', linestyle='None')

    tempdata = trend_data[trend_data['trend'] == 'increasing']
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='^', label='Increasing Trend', linestyle='None')

    tempdata = trend_data[trend_data['trend'] == 'decreasing']
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='v', label='Decreasing Trend', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])
    plt.legend(handels, labels, title='Nitrate Nitrogen Trend', loc='lower right')
    plt.savefig(os.path.join(outdir, 'NO3_N_trends_all_data.png'), dpi=300, bbox_inches='tight')

    plt.close()

    trend_data = pd.DataFrame(index=sites, columns=['trend', 'count', 'northing', 'easting'])
    trend_data.loc[:, 'count'] = data.groupby(site_col).count()
    temp = data.groupby(site_col).aggregate({easting_col: _get_first, northing_col: _get_first})
    trend_data.loc[:, 'lat'] = temp[northing_col]
    trend_data.loc[:, 'lon'] = temp[easting_col]
    for site in sites:
        temp = mann_kendall_obj(data.loc[data[site_col]==site].set_index(datetime_col), data_col=value_col)
        trend_data.loc[site, 'trend'] = temp.trend
    trend_data.to_csv(os.path.join(outdir,'n_trends.csv'))

    # plot overview
    fig, ax = plot_background(gis_dir)
    trend_data = trend_data.loc[trend_data['count'] >= 10]
    tempdata = trend_data[trend_data['trend'] == 'increasing']
    plt.plot(tempdata.lon, tempdata.lat, color='r', marker='^', label='Increasing Trend', linestyle='None')

    tempdata = trend_data[trend_data['trend'] == 'decreasing']
    plt.plot(tempdata.lon, tempdata.lat, color='b', marker='v', label='Decreasing Trend', linestyle='None')

    tempdata = trend_data[trend_data['trend'] == 'no trend']
    plt.plot(tempdata.lon, tempdata.lat, color='grey', marker='o', label='No Trend', linestyle='None')

    plt.xlim([1320000, 1690000])
    plt.ylim([5000000, 5360000])
    plt.tick_params(axis='both', which='both', left='off', labelleft='off', bottom='off', top='off', labelbottom='off')

    handels, labels = ax.get_legend_handles_labels()
    seds = patches.Patch(color='khaki', label='Cenozoic sediments')
    zones = patches.Patch(facecolor='none', edgecolor='k', linewidth=1.5, label='CWMS Zones')
    handels.extend([seds, zones])
    labels.extend([seds.get_label(), zones.get_label()])
    plt.legend(handels, labels, title='Nitrate Nitrogen Trend', loc='lower right')
    plt.savefig(os.path.join(outdir, 'NO3_N_trends_10yr_data.png'), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    outdir = r"C:\Users\MattH\Downloads\test_trends"
    data = pd.read_csv(r"C:\Users\MattH\Downloads\GWQ_A_NO3.csv")
    data.loc[:, 'datetime'] = pd.to_datetime(data['Date'] + ' ' + data['Time'])
    idx = data['Value'].str.contains('<')
    data.loc[:,'Value'] = data.loc[:,'Value'].str.replace('<','').astype(float)
    data.loc[:,'Value'] *= 1/2
    month = [e.month for e in data.datetime]
    idx = np.in1d(month, [9,10,11,12])
    data = data.loc[idx]
    make_n_trends(data, outdir, site_col='Site Name', datetime_col='datetime',
                  value_col='Value', easting_col='Easting', northing_col='Northing')
