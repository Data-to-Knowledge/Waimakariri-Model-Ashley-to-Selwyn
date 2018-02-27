# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 28/08/2017 4:14 PM
"""
from __future__ import division
from pandas import read_csv,read_table
from numpy import float32, int32, sum
base_dir='.'
path = '{}/wel_pkg_base_data.csv'.format(base_dir)
data = read_csv(path, dtype={'layer': int32, 'row': int32, 'col': int32, 'flux': float32})
mult_path = '{}/wel_adj.txt'.format(base_dir) # figure out what the pest control file will populate
multipliers = read_table(mult_path,index_col=0,delim_whitespace=True)

mult_groups = ['pump_c', 'pump_s', 'pump_w', 'sriv', 'n_race', 's_race', 'nbndf']
for group in mult_groups:
    data.loc[data.nsmc_type == group, 'flux'] *= multipliers.loc[group,'value']

add_groups = ['llrzf', 'ulrzf']
for group in add_groups:
    data.loc[data.nsmc_type == group, 'flux'] = multipliers.loc[group,'value']/(data.nsmc_type==group).sum()

g = data.groupby(['layer', 'row', 'col'])
well_ag = g.aggregate({'flux': sum}).reset_index()

outdata = well_ag.loc[:, ['layer', 'row', 'col', 'flux']]
outdata = outdata.rename(columns={'layer': 'k', 'row': 'i', 'col': 'j'})
outdata.loc[:,['k','i','j']] += 1
out_file = r"C:\Users\MattH\Downloads\test_speed.wel" # what will this write as
with open(out_file,'w') as f:
    f.writelines(['11095    740    AUX    IFACE\n',
                  '11095    0    # stress period 0\n'])
outdata.to_csv(out_file,header=None,index=None,sep=' ', mode='a')
