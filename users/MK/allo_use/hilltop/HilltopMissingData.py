# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 15:37:50 2017

@author: MichaelEK
"""

from core.ecan_io.hilltop import rd_hilltop_sites, parse_dsn
from core.ecan_io import write_sql
from os.path import split
from pandas import concat, DataFrame, to_datetime
from configparser import ConfigParser
from os import path, getcwd

from Tkinter import Tk
from tkFileDialog import askdirectory, askopenfilename

Tk().withdraw()

#######################################################
#### Parameters

#ini_file = 'HilltopMissingData.ini'

#### Load in ini parameters
py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))

ini_dir = askopenfilename(initialdir=py_dir, title='Select the ini file', defaultextension='.ini', filetypes=[('ini file', '*.ini')])

ini1 = ConfigParser()
ini1.read([ini_dir])

dsn_path = ini1.get('Input', 'dsn_file')
server = ini1.get('SQLOutput', 'server')
database = ini1.get('SQLOutput', 'database')
sites_table = ini1.get('SQLOutput', 'sites_table')

########################################################
### Run through all hts files

hts_files = parse_dsn(dsn_path)

ht_sites_lst = []
ht_data = DataFrame()
for j in hts_files:
    print(j)
    ## Sites
    sdata = rd_hilltop_sites(j)
    if sdata.empty:
        continue
    base_path, filename = split(j)
    sdata['hts_file'] = filename
    sdata['folder'] = base_path
    ht_sites_lst.append(sdata)

ht_sites = concat(ht_sites_lst).apply(lambda x: x.str.encode('utf-8'), axis=1)
ht_sites.drop('divisor', axis=1, inplace=True)
ht_sites.drop_duplicates(inplace=True)
ht_sites = ht_sites[~ht_sites.mtype.isin(['Regularity'])]

#ht_sites.apply(lambda x: x.str.len(), axis=1).max()
ht_sites.loc[:, 'site'] = ht_sites.loc[:, 'site'].apply(lambda x: x.split('-M')[0])
ht_sites.loc[:, 'start_date'] = to_datetime(ht_sites.loc[:, 'start_date'], errors='coerce')
ht_sites.loc[:, 'end_date'] = to_datetime(ht_sites.loc[:, 'end_date'], errors='coerce')
ht_sites.loc[:, 'folder'] = ht_sites.loc[:, 'folder'].str.replace('\\', '-').str.replace('--', '')
ht_sites = ht_sites.dropna()

## Save to SQl

sites_dtype = {}
for i in ht_sites:
    if i in ['start_date', 'end_date']:
        sites_dtype.update({i: 'DATETIME'})
    else:
        max1 = ht_sites[i].str.len().max()
        sites_dtype.update({i: 'VARCHAR(' + str(max1 + 1) + ')'})

write_sql(server, database, sites_table, ht_sites.copy(), sites_dtype, drop_table=True)














