# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 16:08:44 2017

@author: MichaelEK
"""
from core.allo_use import restr_days, crc_band_flow, ros_freq
from pandas import read_csv
from os import path

############################################
### Parameters

base_dir = r'E:\ecan\shared\projects\hurunui\ros'

sites_csv = 'ros.csv'

ros_csv = 'hurunui_ros_2017-12-06.csv'
months = [1, 2, 3, 4, 10, 11, 12]

###########################################
### Process data

sites = read_csv(path.join(base_dir, sites_csv))

r2 = restr_days(sites.site.tolist(), months=months)

r2.to_csv(path.join(base_dir, ros_csv), index=False)












































