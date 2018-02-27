# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 23/11/2017 11:34 AM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
from core.ecan_io import rd_sql, sql_db
from core.classes.hydro import hydro


def make_n_sheet():
    data = pd.read_excel(
        r"P:\Groundwater\Matt_Hanson\small jobs\N_for_hydrosoc\Copy of Nitrate monitoring data 1995 to 2017.xlsx",
        index_col=0)
    data.loc[:, 'datetime'] = pd.to_datetime(data.loc[:, 'Date'])
    data.loc[:, 'no3'] = [float(str(e).strip('<')) for e in data.loc[:, 'Nitrate Nitrogen']]
    data.loc[:, 'year'] = [e.year for e in data.loc[:, 'datetime']]
    data.loc[:, 'month'] = [e.month for e in data.loc[:, 'datetime']]

    spring2016 = data.loc[(data.year == 2016) & np.in1d(data.month, [9, 10, 11])].reset_index()  # check months
    spring2016 = pd.DataFrame(
        spring2016.groupby('Site_ID').aggregate({'no3': np.mean, 'Site_East': np.mean, 'Site_North': np.mean}))
    spring2017 = data.loc[(data.year == 2017) & np.in1d(data.month, [9, 10, 11])].reset_index()  # check months
    spring2017 = pd.DataFrame(
        spring2017.groupby('Site_ID').aggregate({'no3': np.mean, 'Site_East': np.mean, 'Site_North': np.mean}))

    outdata = pd.merge(spring2016, spring2017, right_index=True, left_index=True, suffixes=[2016, 2017])
    outdata.loc[:, 'ratio'] = outdata.loc[:, 'no32017'] / outdata.loc[:, 'no32016']
    outdata.loc[:, 'rel_change'] = outdata.loc[:, 'ratio'] - 1

    wells = rd_sql(**sql_db.wells_db.well_details).set_index('WELL_NO')
    outdata.loc[:, 'depth'] = [wells.loc[e, 'DEPTH'] for e in outdata.index]
    outdata.to_csv(r"P:\Groundwater\Matt_Hanson\small jobs\N_for_hydrosoc\spring_2017_2016_n.csv")


def make_rainfall():
    h1 = hydro().get_data('precip',sites=r"P:\Groundwater\Matt_Hanson\small jobs\N_for_hydrosoc\canterbury.shp",from_date='2007-01-01').resample('M',fun='sum')
    data = h1.data.loc['atmos_precip_cont_qc'].reset_index()
    data.loc[:, 'month'] = [e.month for e in data.time]
    data.loc[:,'season'] = data.loc[:,'month']
    data.loc[np.in1d(data.season,[6,7,8]),'season'] = 'winter'
    data.loc[:,'year'] = [e.year for e in data.time]
    data = data.groupby(['site','season','year']).aggregate({'data':np.sum}).reset_index()
    data = data.loc[data.season=='winter']
    outdata = data.pivot('site', 'year', 'data')
    outdata.loc[outdata.index,'x'] = [e.x for e in h1.geo_loc.loc[outdata.index, 'geometry']]
    outdata.loc[outdata.index,'y'] = [e.x for e in h1.geo_loc.loc[outdata.index, 'geometry']]
    outdata.to_csv(r"P:\Groundwater\Matt_Hanson\small jobs\N_for_hydrosoc\ecan_rainfall.csv")
    print'done'



if __name__ == '__main__':
    make_rainfall()
    print(';domne')
