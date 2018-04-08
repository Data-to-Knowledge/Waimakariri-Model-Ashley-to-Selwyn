# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/04/2018 8:37 AM
"""

from __future__ import division
from core import env
import pandas as pd


if __name__ == '__main__':
    data = pd.read_excel(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Surface water\Spreadsheets\Waimakariri River N data annual medians.xlsx', sheetname=2)
    data.loc[:,'year'] =pd.to_datetime(data.Date, errors='coerce').dt.year
    outdata = data.groupby(['Site','year']).describe().loc[:,'Nitrate/nitrite (g/m3 N)']
    outdata.to_csv(r'\\gisdata\projects\SCI\Groundwater\Waimakariri\Surface water\Spreadsheets\waimak_annual_medians_niwa.csv')

