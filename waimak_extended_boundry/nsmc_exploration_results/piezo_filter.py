# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 12:40 PM
"""

from __future__ import division
import pandas as pd
from waimak_extended_boundry.supporting_data_analysis.all_well_layer_col_row import get_all_well_row_col

def extract_peizo_filter_raw_data(outpath):
    piezo_data = pd.read_csv(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\2017 piezo survey\Survey results\DataforNSMC_Piezo.csv",index_col=1)
    all_wells = get_all_well_row_col()
    outdata = pd.merge(pd.DataFrame(piezo_data.loc[:,'Sept17PiezoElev']),
                       all_wells.loc[:,['layer', 'nztmx', 'nztmy']],
                       left_index=True,right_index=True)
    outdata.rename({'Sept17PiezoElev': 'obs'},inplace=True)
    outdata.to_csv(outpath)

# this was done by brioch, this file is depreciated

if __name__ == '__main__':
    extract_peizo_filter_raw_data(r"C:\Users\MattH\Desktop\piezo_data_for_filter.csv")
    print('done')