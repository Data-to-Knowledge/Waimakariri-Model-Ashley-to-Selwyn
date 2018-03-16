# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 16/03/2018 4:24 PM
"""

from __future__ import division
from core import env
import pandas as pd
import os

if __name__ == '__main__':
    sites = [
        'ashely_gorge',
        'ashely_sh1',
        'whole_estuary',
        'saltwater_end',
        'waikuku_end',
        'taranaki_end',
    ]
    landuse_types = [u'',
                     u'Arable',
                     u'DairyFarm',
                     u'DairySupport',
                     u'Forest-Tussock',
                     u'Horticulture',
                     u'Lifestyle',
                     u'NotFarm',
                     u'Other-inclGolf',
                     u'Pigs',
                     u'Sheep-Beef-Deer',
                     u'SheepBeef-Hill',
                     u'Unknown']

    base_dir = r"C:\Users\MattH\Downloads\test_nload_stuffs"

    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_landuse.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_total_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_landuse.csv'))

    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_landuse.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'area_in_catch']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_sizes_landuse.csv'))

    landuse_types = ['0-10 ha',
                     '10-20 ha',
                     '20-50 ha',
                     '>50 ha', ]
    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_size.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_total_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_size.csv'))

    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_size.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'area_in_catch']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_sizes_size.csv'))

    #todo  winter vs irrigation
    landuse_types = [u'',
                     u'Arable',
                     u'DairyFarm',
                     u'DairySupport',
                     u'Forest-Tussock',
                     u'Horticulture',
                     u'Lifestyle',
                     u'NotFarm',
                     u'Other-inclGolf',
                     u'Pigs',
                     u'Sheep-Beef-Deer',
                     u'SheepBeef-Hill',
                     u'Unknown']
    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_landuse.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_forage_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_winter_landuse.csv'))

    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_landuse.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_irrigation_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_irrig_landuse.csv'))

    landuse_types = ['0-10 ha',
                     '10-20 ha',
                     '20-50 ha',
                     '>50 ha', ]
    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_size.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_forage_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_winter_size.csv'))

    outdata_sites = pd.DataFrame(columns=sites, index=landuse_types)
    for site in sites:
        tempdata = pd.read_csv(os.path.join(base_dir, site + '_grouped_size.csv'), index_col=0)
        for lt in landuse_types:
            try:
                outdata_sites.loc[lt, site] = tempdata.loc[lt, 'new_PA_irrigation_load_prorata']
            except KeyError:
                pass
    outdata_sites.to_csv(os.path.join(base_dir, 'ashleytribs_overview_paloads_irrig_size.csv'))
