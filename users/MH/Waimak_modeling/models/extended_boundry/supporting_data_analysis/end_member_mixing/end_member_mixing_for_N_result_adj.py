# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/05/2018 11:26 AM
"""

from __future__ import division
from core import env
import pandas as pd
from waimak_4_endmember_mixing import mc_calc_end_members

if __name__ == '__main__':
    end_mean_o18 = {'inland': -8.76, 'coastal': -8.00, 'river': -9.25, 'eyre': -8.90}
    end_sd_o18 = {'inland': 0.21, 'coastal': 0.64, 'river': 0.22, 'eyre': 0.54}

    end_mean_cl = {'inland': 10.38, 'coastal': 25.43, 'river': 1.05, 'eyre': 3.67}
    end_sd_cl = {'inland': 4.21, 'coastal': 9.49, 'river': 0.15, 'eyre': 0.29}

    targets = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\emma_for_n_adjustment\cl_o18_data_grouped.csv",
        index_col=0)
    sites = {}
    for site in targets.index:
        sites[site] = {'o18_lower': targets.loc[site, 'mean_o18'] - targets.loc[site, 'std_o18'],
                       'o18_upper': targets.loc[site, 'mean_o18'] + targets.loc[site, 'std_o18'],
                       'cl_lower': targets.loc[site, 'mean_cl'] - targets.loc[site, 'std_cl'],
                       'cl_upper': targets.loc[site, 'mean_cl'] + targets.loc[site, 'std_cl']}
    mc_calc_end_members(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\emma_for_n_adjustment\4_endmembers",
        sites, end_mean_o18, end_sd_o18, end_mean_cl, end_sd_cl,n=1000, export_raw=True)
    print('done')
