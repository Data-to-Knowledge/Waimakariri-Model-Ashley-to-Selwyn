# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/01/2018 9:01 AM
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
import os


def load_sd(path, sd_version):
    out = pd.read_csv(os.path.join(path, 'NsmcBase_extract_{}.csv'.format(sd_version)), skiprows=1, index_col=0)
    out = out.loc[:, [u'custmaindrain_swaz', u'cust_swaz', u'eyre_swaz',
                      u'n7drain_swaz', u'cam_swaz', u'courtenay_swaz', u'greigs_swaz',
                      u'kaiapoi_swaz', u'kairaki_swaz', u'northbrook_swaz', u'ohoka_swaz',
                      u'saltwater_swaz', u'southbrook_swaz', u'taranaki_swaz',
                      u'waikuku_swaz', u'ashley_swaz', u'waimakupper_swaz',
                      u'waimaklower_swaz']]
    return out


if __name__ == '__main__':
    high_base =r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\high_s\NsmcBase_2017_12_22_data"
    med_base =r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\med_s\NsmcBase_2017_12_26_data"
    low_base =r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\low_s\NsmcBase_2017_12_26_data"
    high = {}
    med = {}
    low = {}

    for sdv in ['sd7', 'sd30', 'sd150']:
        high[sdv] = load_sd(high_base, sdv)
        med[sdv] = load_sd(med_base, sdv)
        low[sdv] = load_sd(low_base, sdv)

        lows = low[sdv].values.flatten()
        highs = high[sdv].values.flatten()
        meds = med[sdv].values.flatten()

        fig, ax = plt.subplots()
        ax.scatter(lows, highs)
        ax.plot(ax.get_ylim(),ax.get_ylim(),ls = '--')
        ax.set_ylabel('high')
        ax.set_xlabel('low')
        ax.set_title('{} low vs high'.format(sdv))

        fig, ax = plt.subplots()
        ax.scatter(meds, highs)
        ax.plot(ax.get_ylim(),ax.get_ylim(),ls = '--')
        ax.set_ylabel('high')
        ax.set_xlabel('med')
        ax.set_title('{} med vs high'.format(sdv))

        fig, ax = plt.subplots()
        ax.scatter(lows, meds)
        ax.plot(ax.get_ylim(),ax.get_ylim(),ls = '--')
        ax.set_ylabel('med')
        ax.set_xlabel('low')
        ax.set_title('{} low vs med'.format(sdv))
    plt.show()
