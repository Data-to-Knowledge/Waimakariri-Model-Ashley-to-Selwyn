# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 24/08/2017 2:42 PM
"""

from __future__ import division
from waimak_extended_boundry import _get_wel_spd_v1
import pandas as pd
import numpy as np


if __name__ == '__main__':
    well_data = _get_wel_spd_v1()
    well_data = well_data.loc[well_data.type=='well']

    mike = pd.read_hdf(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\sd_est_all_mon_vol.h5")
    mike = mike.loc[(mike.time >= pd.datetime(2008, 1, 1)) & (mike.take_type == 'Take Groundwater')]
    mike.loc[:, 'd_in_m'] = mike.time.dt.daysinmonth
    data = mike.groupby('wap').aggregate({'mon_allo_m3': np.sum, 'usage_est': np.sum,'crc': ','.join, 'd_in_m': np.sum})
    data.loc[:, 'flux'] = data.loc[:, 'mon_allo_m3'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days

    # assume that at the outside option any well can be at 0 or CAV (and that this range represents the 99.7 possibliity)

    well_data.loc[:,'consented_flux'] = data.loc[well_data.index,'mon_allo_m3']
    well_data.loc[:,'est_flux'] = data.loc[well_data.index,'usage_est']

    for cwms in set(well_data.cwms):
        if pd.isnull(cwms):
            continue
        print('')
        print(cwms)
        temp = well_data.loc[well_data.cwms==cwms]
        print(((temp.flux)**2).sum()**0.5/(temp.flux.sum()*-1))
        per_use = temp.est_flux.sum()/temp.consented_flux.sum()
        print('percentage use: {}'.format(per_use))

