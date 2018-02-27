# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 09:12:25 2017

@author: MichaelEK
"""

from core.ecan_io.hilltop import ht_sites, ht_get_data
from pandas import concat, merge

###############################################
### Parameters

hts1 = r'H:\Data\Archive\Boraman2016-17.hts'
hts2 = r'H:\Data\Telemetry\Boraman.hts'

sites_out = r'E:\ecan\local\Projects\requests\Ilja\2017-09-18\boraman_sites_comparison.csv'

##############################################


sites1 = ht_sites(hts1)
sites2 = ht_sites(hts2)

sites_all = merge(sites1, sites2, on='site', how='outer', suffixes=('_archive', '_tel'))

sites_all.to_csv(sites_out, index=False, header=True)



























