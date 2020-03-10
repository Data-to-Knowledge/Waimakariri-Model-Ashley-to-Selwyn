# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/06/2018 12:02 PM
"""

from __future__ import division
from glob import glob
import pandas as pd
import os


if __name__ == '__main__':
    well_num = 'BW24/0074'
    paths = glob(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\*\waimakariri_zone\raw_model_data\*\stocastic_set_wells.csv")
    names = [os.path.basename(os.path.dirname(e)) for e in paths]
    outdata = pd.DataFrame(index=names, columns=[u'count', u'mean', u'std', u'min', u'5%', u'25%', u'50%', u'75%',
       u'95%', u'max', u'Zone', u'Zone_1', u'zone_2', u'private_public'])
    for n, p in zip(names, paths):
        temp = pd.read_csv(p,index_col=0)
        outdata.loc[n] = temp.loc[well_num]
    outdata.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\Tuahiwi_Marae_supply_well.csv")