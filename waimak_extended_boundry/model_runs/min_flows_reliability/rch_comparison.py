# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/10/2017 5:14 PM
"""

from __future__ import division
import numpy as np
from waimak_extended_boundry.model_and_NSMC_build.m_packages.rch_packages import _get_rch
from waimak_extended_boundry.model_run_tools import get_lsrm_base_array, get_zone_array_index
from copy import deepcopy
import itertools
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import os

def get_rch_through_time():
    """
    exploritory plots
    :return:
    """
    periods = range(2010, 2100, 20)
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    senarios = ['pc5', 'nat', 'current']

    w_idx = get_zone_array_index(['waimak'])
    outdata = pd.DataFrame(list(itertools.product(periods, rcms, rcps, amalg_types,senarios)),columns=['period','rcm','rcp','at','senario'])
    for i,per, rcm, rcp, at, sen in outdata.itertuples(name=None):
        new_rch = get_lsrm_base_array(sen,rcp,rcm,per,at)
        outdata.loc[i,'rch'] = np.nansum(new_rch[w_idx])*200*200/86400
    temp = pd.DataFrame(list(itertools.product([1980], rcms, ['RCPpast'], amalg_types, senarios)),columns=['period','rcm','rcp','at','senario'])
    for i,per, rcm, rcp, at, sen in temp.itertuples(name=None):
        new_rch = get_lsrm_base_array(sen, rcp, rcm, per, at)
        temp.loc[i,'rch'] = np.nansum(new_rch[w_idx])*200*200/86400
    temp2 = deepcopy(temp)
    temp.loc[:,'period'] = 1990
    temp2.loc[:,'period'] = 1990
    temp.loc[:,'rcp'] = 'RCP4.5'
    temp2.loc[:,'rcp'] = 'RCP8.5'
    outdata = pd.concat((outdata,temp,temp2))

    temp = outdata.loc[(outdata['at'] == 'period_mean')]
    g = sns.FacetGrid(temp,row='senario',col='rcp')
    g.map_dataframe(sns.pointplot,x='period',y='rch',hue='rcm')
    plt.show()
    print'done'

def model_average_rcp_past():
    """
    didn't finish this
    :return:
    """
    paths = glob(r"K:\niwa_netcdf\lsrm\lsrm_results\water_year_means\wym_RCPpast_*_80perc.h5")
    outdata = {}
    for path in paths:
        rcm = os.path.basename(path).split('_')[-2]
        temp = pd.read_hdf(path)
        temp.loc[:,'year'] = [e.year for e in temp.time]
        g = temp.groupby('year').aggregate({'total_drainage': np.sum})

    raise NotImplementedError




if __name__ == '__main__':
    get_rch_through_time()

    w_idx = get_zone_array_index('waimak')
    old_rch = _get_rch()
    new_rch = get_lsrm_base_array('current',None,None,None,'mean')
    new_rch[~np.isfinite(new_rch)] = 0

    print 'done'