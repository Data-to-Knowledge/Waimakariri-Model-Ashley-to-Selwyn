# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 5/07/2017 3:55 PM
"""

from __future__ import division
from core import env
from core.ecan_io import rd_sql, sql_db
import pandas as pd
import rasterio


if __name__ == '__main__':


    with open(r"C:\Users\MattH\Downloads\old_targets.txt") as f:
        wells = f.readlines()
    wells = [e.strip() for e in wells]
    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details_org = well_details_org.set_index('WELL_NO')

    rb = rasterio.open(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\tops.tif")

    outdata = pd.DataFrame(index=wells)
    for well in outdata.index:
        x = outdata.loc[well,'x'] = well_details_org.loc[well,'NZTMX']
        y = outdata.loc[well,'y'] = well_details_org.loc[well,'NZTMY']
        outdata.loc[well,'depth'] = well_details_org.loc[well,'DEPTH']
        outdata.loc[well,'dem_ground'] = list(rb.sample([(x, y)]))[0][0]
        ref = outdata.loc[well,'ref_level'] = well_details_org.loc[well,'REFERENCE_RL']
        ref_ground = outdata.loc[well,'ref_ground'] = well_details_org.loc[well,'GROUND_RL']
        if pd.isnull(ref_ground):
            ref_ground = 0
        outdata.loc[well, 'ground_level'] = ref + ref_ground
    outdata.loc[:,'diff'] = outdata.loc[:,'ground_level'] - outdata.loc[:,'dem_ground']
    outdata.to_csv(r"C:\Users\MattH\Downloads\old_targets_ground_comparison2.txt")

    print('done')




