# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 16:48:08 2017

@author: KateSt
"""
from core import env
from core.ecan_io import rd_sql, sql_db
from pandas import read_csv, merge

test = rd_sql(server='SQL2012PROD03', database='WUS', table='vConsent_Usage_2015')

consents_69635 = read_csv(filepath_or_buffer = r'S:\Surface Water\shared\projects\otop\naturalisation\69635 - Consents.csv')

consents_69635_1 = consents_69635.drop_duplicates(subset='crc')

consents_69635_2 = consents_69635_1['crc']


consents_69635_usage = test.isin(consents_69635_2)

sd1.to_csv(path_or_buf=r'S:\Surface Water\shared\projects\otop\naturalisation\sd1.csv')

sd1a.to_csv(path_or_buf=r'S:\Surface Water\shared\projects\otop\naturalisation\sd1a.csv')