# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/08/2017 3:15 PM
"""

from __future__ import division
import numpy as np
import pandas as pd
from core.ecan_io import rd_sql, sql_db

mike = pd.read_hdf(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\sd_est_all_mon_vol.h5")
mike = mike.loc[(mike.time >= pd.datetime(2008, 1, 1)) & (mike.take_type == 'Take Groundwater')]
mike.loc[:, 'd_in_m'] = mike.time.dt.daysinmonth
data = mike.groupby('wap').aggregate({'usage_est': np.sum, 'mon_allo_m3': np.sum, 'crc': ','.join, 'd_in_m': np.sum})
data.loc[:, 'flux'] = data.loc[:, 'usage_est'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days
data.loc[:, 'flux_cav'] = data.loc[:, 'mon_allo_m3'] / (mike.time.max() - pd.datetime(2007, 12, 31)).days

well_details = rd_sql(**sql_db.wells_db.well_details)
well_details = well_details.set_index('WELL_NO')
out_data = pd.merge(data, pd.DataFrame(well_details.loc[:, 'WMCRZone']), left_index=True, right_index=True)
out_data = out_data.loc[np.in1d(out_data.WMCRZone, [7, 8])]

temp = out_data.flux/out_data.flux_cav
temp2 = temp[temp<=1]

print 'done'