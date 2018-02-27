"""
coding=utf-8
Author: matth
Date Created: 29/06/2017 8:51 AM

Fixed!!

"""

from __future__ import division
from core import env
from core.ecan_io import rd_sql, sql_db
import pandas as pd
from core.classes.hydro import hydro
import rasterio


def get_mean_water_level():
    well_details_org = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details_org[(well_details_org['WMCRZone'] == 4) | (well_details_org['WMCRZone'] == 7) |
                                (well_details_org['WMCRZone'] == 8)]  # keep only waimak selwyn and chch zones
    well_details = well_details[well_details['Well_Status'] == 'AE']
    well_details = well_details[pd.notnull(well_details['DEPTH'])]

    well_details.loc[:,'WELL_NO'] = [e.strip() for e in well_details.loc[:,'WELL_NO']]
    well_details = well_details.set_index('WELL_NO')
    # why are missing sites coming up ask mike?

    data = hydro().get_data(mtypes=['gwl_m'], sites=list(well_details.index))
    temp = data.groupby(level=['mtype', 'site']).describe()[['min', '25%', '50%', '75%', 'mean', 'max', 'count']].round(2)
    sites = list(temp.loc['gwl_m'].index[temp.loc['gwl_m']['count'] >=5])

    data = data.loc['gwl_m',sites].reset_index()
    data = data.set_index('site')
    data['month'] = [e.month for e in data.time]
    data['year'] = [e.year for e in data.time]


    out_data = pd.DataFrame(index=set(data.index))

    for well in out_data.index:
        out_data.loc[well,'nztmx'] = well_details.loc[well,'NZTMX']
        out_data.loc[well,'nztmy'] = well_details.loc[well,'NZTMY']
        out_data.loc[well,'depth'] = well_details.loc[well,'DEPTH']


gwl_m_dict = {'server': prod_server05, 'db': 'Wells', 'table': 'DTW_READINGS', 'site_col': 'WELL_NO', 'time_col': 'DATE_READ', 'data_col': 'DEPTH_TO_WATER', 'qual_col': None}

t1 = rd_sql(gwl_m_dict['server'], gwl_m_dict['db'], gwl_m_dict['table'], ['WELL_NO', 'DATE_READ', 'DEPTH_TO_WATER'], where_col='WELL_NO', where_val=['L36/0124'])
t1['WELL_NO'].unique()



t2 = [u'L36/0124 ', u'M35/0132 ', u'M35/0698 ', u'M35/1003 ', u'M35/1103 ', u'M35/4682 ', u'M35/5918 ', u'M35/6656 ', u'M36/1318 ', u'M36/1421 ', u'M36/4050 ', u'M36/4227 ', u'M36/5248 ']

sites2 = [str(i) for i in list(well_details.index)]

data = hydro().get_data(mtypes=['gwl_m'], sites=sites2)













