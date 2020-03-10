# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/08/2017 1:41 PM
"""

from __future__ import division
import numpy as np
import pandas as pd
from core.ecan_io import rd_sql, sql_db
from waimak_extended_boundry import smt
from waimak_extended_boundry import get_wel_spd

well_data = get_wel_spd(0)
well_data = well_data.loc[(well_data.type=='well')]

mike = pd.read_hdf(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\sd_est_all_mon_vol.h5")
mike = mike.loc[(mike.time>=pd.datetime(2008,1,1)) & (mike.take_type == 'Take Groundwater')]
mike.loc[:,'d_in_m'] = mike.time.dt.daysinmonth
data = mike.groupby('wap').aggregate({'usage_est': np.sum, 'crc': ','.join,'d_in_m':np.sum})
data.loc[:,'flux'] = data.loc[:,'usage_est']/(mike.time.max() - pd.datetime(2007,12,31)).days

well_details = rd_sql(**sql_db.wells_db.well_details)
well_details = well_details.set_index('WELL_NO')
s_wai_wells = pd.merge(data,pd.DataFrame(well_details.loc[:,'WMCRZone']), left_index=True, right_index=True)
s_wai_wells = s_wai_wells.loc[np.in1d(s_wai_wells.WMCRZone,[7,8,4])]
s_wai_wells.loc[:,'flux'] *= -1

temp = smt.get_well_postions(np.array(s_wai_wells.index), one_val_per_well=True, raise_exct=False)
s_wai_wells['layer'], s_wai_wells['row'], s_wai_wells['col'] = temp
no_flow = smt.get_no_flow()
for i in s_wai_wells.index:
    layer, row, col = s_wai_wells.loc[i, ['layer', 'row', 'col']]
    if any(pd.isnull([layer, row, col])):
        continue
    if no_flow[layer, row, col] == 0:  # get rid of non-active wells
        s_wai_wells.loc[i, 'layer'] = np.nan

s_wai_wells = s_wai_wells.dropna(subset=['layer', 'row', 'col'])

flux_stats = pd.DataFrame(index=['Ashley sub-zone', 'Cust sub-zone', 'Eyre sub-zone','Waimakariri CWMS',
                                 'CHCH-WM CWMS','Selwyn CWMS'],
                          columns=['model','mike'])
arrays = {}
for name, dat in zip(['model','mike'],[well_data,s_wai_wells]):
    # waimak = 4, chch_wm = 7, selwyn=8
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    # 1 = ashley, 2 = eyre, 3 = cust
    sub_zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/Sub_cwms_zones.shp".format(smt.sdp),'id_num')
    flux_array = smt.df_to_array(dat,'flux')
    arrays[name] = flux_array
    for zone_id, zone_name in zip([1, 3, 2,
                                   4, 7, 8],
                                  ['Ashley sub-zone', 'Cust sub-zone', 'Eyre sub-zone',
                                   'Waimakariri CWMS', 'CHCH-WM CWMS','Selwyn CWMS']):
        if zone_id > 3:
            id_matrix = zones
        else:
            id_matrix = sub_zones

        flux_stats.loc[zone_name,name] = np.nansum(flux_array[np.isclose(id_matrix,zone_id)])
flux_stats = flux_stats/86400

print(flux_stats)

print 'done'