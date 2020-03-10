# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 2:21 PM
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
from waimak_extended_boundry import get_all_well_row_col

if __name__ == '__main__':
    highs = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\high_s\NsmcBase_2017_12_22_data\NsmcBase_extract_sd150.csv", index_col=0, skiprows=1)
    meds = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\med_s\NsmcBase_2017_12_26_data\NsmcBase_extract_sd150.csv", index_col=0, skiprows=1)
    lows = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\low_s\NsmcBase_2017_12_26_data\NsmcBase_extract_sd150.csv",index_col=0, skiprows=1)
    theis = pd.read_excel(r"\\gisdata\Projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\stream_depletion\NCS_NumericalvsAnalyticalModelComparison.xlsx", sheetname='Theis Matt Smith 180118', index_col=0)
    all_wells = get_all_well_row_col()
    sumkeys =[
        'custmaindrain_swaz',
        'cust_swaz',
        'eyre_swaz',
        'n7drain_swaz',
        'cam_swaz',
        'courtenay_swaz',
        'greigs_swaz',
        'kaiapoi_swaz',
        'kairaki_swaz',
        'northbrook_swaz',
        'ohoka_swaz',
        'saltwater_swaz',
        'southbrook_swaz',
        'taranaki_swaz',
        'waikuku_swaz',
        'ashley_swaz',

    ]

    outdata = pd.DataFrame(index=highs.index,columns=['high_s','med_s','low_s'])
    outdata.loc[:,'high_s'] = highs.loc[:,sumkeys].sum(axis=1)
    outdata.loc[:,'med_s'] = meds.loc[:,sumkeys].sum(axis=1)
    outdata.loc[:,'low_s'] = lows.loc[:,sumkeys].sum(axis=1)

    outdata = pd.merge(outdata, pd.DataFrame(theis.loc[:,'SD1_150']), right_index=True, left_index=True)
    outdata = pd.merge(outdata, all_wells, right_index=True, left_index=True)
    outdata = outdata.sort_values(by='depth')
    outdata.to_csv(r"C:\Users\MattH\Downloads\test_sd_comp.csv")

    fig,ax = plt.subplots()

    ax.scatter(outdata.depth, outdata.med_s, label='med_s', c='g')
    ax.scatter(outdata.depth, outdata.high_s, label='high_s', c='orange')
    ax.scatter(outdata.depth, outdata.low_s, label='low_s', c='b')
    ax.scatter(outdata.depth, outdata.SD1_150, label='thies', c='r')
    ax.legend()

    fig,ax = plt.subplots()

    ax.scatter(outdata.SD1_150, outdata.low_s, label='xmatts, ylows', c='g')
    ax.legend()

    fig,ax = plt.subplots()
    x = range(len(outdata.index))
    ax.scatter(x, outdata.med_s, label='med_s', c='g')
    ax.scatter(x, outdata.high_s, label='high_s', c='orange')
    ax.scatter(x, outdata.low_s, label='low_s', c='b')
    ax.scatter(x, outdata.SD1_150, label='thies', c='r')
    ax.legend()

    fig,ax = plt.subplots()
    x = range(len(outdata.index))
    ax.scatter(x, outdata.med_s, label='med_s', c='b')
    ax.scatter(x, outdata.SD1_150, label='thies', c='r')
    ax.legend()

    fig,ax = plt.subplots()
    x = range(len(outdata.index))
    ax.scatter(x, outdata.med_s, label='med_s', c='g')
    ax.scatter(x, outdata.high_s, label='high_s', c='orange')
    ax.scatter(x, outdata.low_s, label='low_s', c='b')
    ax.legend()
    plt.show()
    #  talk with zeb about this... we bailed