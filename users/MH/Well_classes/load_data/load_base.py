"""
Author: matth
Date Created: 28/02/2017 1:37 PM
"""

from __future__ import division
from core.ecan_io import rd_sql, sql_db

# (MH) shift to Well details on AGWQL  this involves mapping the CWMS zones codes to names
def get_base_well_data_df(well_list):
    alldata = rd_sql(**sql_db.wells_db.base)
    out_data = alldata[alldata.WELL_NO.isin(well_list)]

    return out_data