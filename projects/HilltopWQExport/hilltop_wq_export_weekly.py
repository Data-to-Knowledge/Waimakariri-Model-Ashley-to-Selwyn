# -*- coding: utf-8 -*-
"""
Created on Tue Oct 03 08:51:46 2017

@author: MichaelEK
"""
######################################
### Change pythonpath and add init

import sys
from os import path, getcwd

base_py_path = r'D:\Executables\PythonScripts'
proj_name = 'HilltopWQExport'
py_path = path.join(base_py_path, proj_name)

sys.path.append(py_path)

init1 = path.join(py_path, '__init__.py')
fh = open(init1,'a+')

#####################################
### Load Libraries

from time import sleep, time
from numpy import ceil, array_split
from pandas import concat, DataFrame, merge, to_datetime
from configparser import ConfigParser
from datetime import datetime
from sqlalchemy.types import String
from core.ecan_io.hilltop import rd_ht_wq_data, rd_hilltop_sites
from core.ecan_io.mssql import to_mssql, create_mssql_table, del_mssql_table_rows, rd_sql

start_run = time()

#############################################
### Parameters
## Import
print('load parameters')

py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))

#py_dir = r'E:\ecan\git\HydroPandas\users\MK\tsdata_exports'

ini1 = ConfigParser()
ini1.read([path.join(py_dir, path.splitext(__file__)[0] + '.ini')])

#ini1.read([path.join(py_dir, 'hilltop_wq_export_weekly' + '.ini')])

hts_files = [i.strip() for i in ini1.get('Input', 'hts_files').split(',')]

sample_params_list = ['Project', 'Cost Code', 'Technician', 'Sample ID', 'Sample Comment', 'Field Comment', 'Sample Appearance', 'Sample Colour', 'Sample Odour', 'Water Colour', 'Water Clarity']
mtype_params_list = ['Lab Method', 'Lab Name']

sample_params_list = [i.strip() for i in ini1.get('Input', 'sample_params_list').split(',')]
mtype_params_list = [i.strip() for i in ini1.get('Input', 'mtype_params_list').split(',')]

#sites_chunk = 100
sites_chunk = int(ini1.get('Input', 'sites_chunk'))

## Export

sites_dtype = {'SiteID': 'VARCHAR(19)', 'MeasurementType': 'VARCHAR(59)', 'Units': 'VARCHAR(29)', 'FromDate': 'DATETIME', 'ToDate': 'DATETIME'}

samples_dtype = {'SiteID': 'VARCHAR(19)', 'CollectionTime': 'DATETIME', 'Param': 'VARCHAR(59)', 'Value': 'VARCHAR(255)'}

mtypes_dtype = {'SiteID': 'VARCHAR(19)', 'MeasurementType': 'VARCHAR(59)', 'CollectionTime': 'DATETIME', 'Param': 'VARCHAR(59)', 'Value': 'VARCHAR(255)'}

col_name_convert = {'site': 'SiteID', 'mtype': 'MeasurementType', 'unit': 'Units', 'start_date': 'FromDate', 'end_date': 'ToDate', 'time': 'CollectionTime'}
col_name_convert2 = {col_name_convert[i]: i for i in col_name_convert}


server = str(ini1.get('Output', 'server'))
database = str(ini1.get('Output', 'database'))
sites_table = str(ini1.get('Output', 'sites_table'))
mtypes_table = str(ini1.get('Output', 'mtypes_table'))
samples_table = str(ini1.get('Output', 'samples_table'))

sites_pkeys = ['SiteID', 'MeasurementType']
samples_pkeys = ['SiteID', 'Param', 'CollectionTime']
mtypes_pkeys = ['SiteID', 'MeasurementType', 'Param', 'CollectionTime']

#sql_log = {'server': 'SQL2014DEV01', 'database': 'Apps_Trace_Log', 'table': 'Log'}
log_server = str(ini1.get('Output', 'log_server'))
log_database = str(ini1.get('Output', 'log_database'))
log_table = str(ini1.get('Output', 'log_table'))

### Start processing

try:

    ###############################################
    ### Extract data

    mtype_list = list(mtype_params_list)
    mtype_list.append('data')

    ## Extract site data
    print('Extract site data')
    sites_dict = {}
    for hts in hts_files:
        sites_info1 = rd_hilltop_sites(hts)
        sites_info2 = sites_info1.drop(['data_source', 'divisor'], axis=1)
        sites_info2.rename(columns=col_name_convert, inplace=True)
        sites_dict.update({hts: sites_info2})

    sites_all = concat(list(sites_dict.values())).reset_index(drop=True)
    sites_all['FromDate'] = to_datetime(sites_all['FromDate'])
    sites_all['ToDate'] = to_datetime(sites_all['ToDate'])
    sites_all['MeasurementType'] = sites_all['MeasurementType'].str.replace("'", '')
    sites_all['Units'] = sites_all['Units'].str.replace("'", '')

    # Check for duplicate sites/mytpes
    dups = sites_all[sites_all.duplicated(sites_pkeys)]

    if not dups.empty:
        raise ValueError('There are duplicate sites and measurement types')

    ### Determine if tables exist and create if needed
    print('Determine if tables exist and create if needed')

    today1 = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    sites_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + sites_table + "', 'U')").loc[0][0] is None
    samples_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + samples_table + "', 'U')").loc[0][0] is None
    mtypes_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + mtypes_table + "', 'U')").loc[0][0] is None

    if samples_bool1:
        create_mssql_table(server, database, samples_table, dtype_dict=samples_dtype, primary_keys=samples_pkeys)

    if mtypes_bool1:
        create_mssql_table(server, database, mtypes_table, dtype_dict=mtypes_dtype, primary_keys=mtypes_pkeys)

    if sites_bool1:
        create_mssql_table(server, database, sites_table, dtype_dict=sites_dtype, primary_keys=sites_pkeys)

    ### Upload sites table to MSSQL
    print('Compare with existing sites table and upload to MSSQL if necessary')
    ## Check if data is the same - sites
    old_site_data = rd_sql(server, database, sites_table)
    old_site_data['FromDate'] = to_datetime(old_site_data['FromDate'])
    old_site_data['ToDate'] = to_datetime(old_site_data['ToDate'])

    site_combo = merge(sites_all, old_site_data, on=sites_pkeys, how='outer', suffixes=['_1', '_2'], indicator=True)
    bad_old_site = site_combo[site_combo._merge == 'right_only']
    if not bad_old_site.empty:
        del_mssql_table_rows(server, database, sites_table, pk_df=bad_old_site[['SiteID', 'MeasurementType']])
        site_combo = site_combo[site_combo._merge != 'right_only']
        print('These sites were removed')
        print(bad_old_site)

    site_combo.set_index(sites_pkeys, inplace=True)
    new_set = site_combo[['Units_1', 'FromDate_1', 'ToDate_1']].copy()
    new_set.columns = ['Units', 'FromDate', 'ToDate']
    old_set = site_combo[['Units_2', 'FromDate_2', 'ToDate_2']].copy()
    old_set.columns = ['Units', 'FromDate', 'ToDate']

    wq_data_site1 = new_set.loc[(new_set[['FromDate', 'ToDate']] != old_set[['FromDate', 'ToDate']]).any(1)].reset_index().copy()
    if not wq_data_site1.empty:
        if not sites_bool1:
            del_mssql_table_rows(server, database, sites_table, pk_df=wq_data_site1[['SiteID', 'MeasurementType']])
        to_mssql(wq_data_site1, server, database, sites_table)
        print(str(len(wq_data_site1)) + ' sites updated')
        print(wq_data_site1)

    ### Extract ts data
    print('Extract ts data')
    for hts in hts_files:
        print(hts)
        hts_sites_all = sites_dict[hts].copy()

        sites1 = hts_sites_all.SiteID.unique()
        n_chunks = ceil(len(sites1) / float(sites_chunk))
        sites2 = array_split(sites1, n_chunks)

        ## Chunk out the data and into MSSQL
        print('There are ' + str(n_chunks) + ' chunks of sites to extract')
        for hts_sites in sites2:

            ## Pull out data and reorganize - due to the Hilltop COM being buggy, the call needs to be tried multiple times
            timer = 5
            while timer > 0:
                try:
                    wq_data = rd_ht_wq_data(hts, sites=hts_sites.tolist(), output_site_data=False, sample_params=sample_params_list, mtype_params=mtype_params_list)
                    break
                except Exception as err:
                    err1 = err
                    timer = timer - 1
                    if timer == 0:
                        raise ValueError(err1)
                    else:
                        print(err1)
                        sleep(3)

            wq_data1 = wq_data.melt(['site', 'mtype', 'time'], var_name='Param', value_name='Value')
            wq_data1['Value'] = wq_data1['Value'].astype(str)

            wq_data_mtype = wq_data1[wq_data1['Param'].isin(mtype_list)].copy()

            wq_data_sample = wq_data1[wq_data1['Param'].isin(sample_params_list)]
            wq_data_sample = wq_data_sample.drop('mtype', axis=1).drop_duplicates(['site', 'Param', 'time'])
            wq_data_sample['Value'] = wq_data_sample['Value'].str.strip()
            wq_data_sample = wq_data_sample[~((wq_data_sample['Value'] == '') | (wq_data_sample['Value'] == '..'))]

            # Rename columns
            wq_data_sample.rename(columns=col_name_convert, inplace=True)
            wq_data_mtype.rename(columns=col_name_convert, inplace=True)

            # Corrections to mtypes
            wq_data_mtype['MeasurementType'] = wq_data_mtype['MeasurementType'].str.replace("'", '')

            ## SQL upload
            # Check if data is the same - samples
            old_sample_data = rd_sql(server, database, samples_table, where_col={'SiteID': hts_sites.tolist()})
            old_sample_data.rename(columns={'Value': 'old'}, inplace=True)
            sample_combo = merge(wq_data_sample, old_sample_data, on=samples_pkeys, how='outer')
            bad_old_sample = sample_combo[sample_combo.Value.isnull()]
            if not bad_old_sample.empty:
                del_mssql_table_rows(server, database, samples_table, pk_df=bad_old_sample[['SiteID', 'Param', 'CollectionTime']])
                sample_combo = sample_combo[sample_combo.Value.notnull()]

            wq_data_sample1 = sample_combo[sample_combo.Value != sample_combo.old].drop('old', axis=1).copy()

            if not wq_data_sample1.empty:
                if not samples_bool1:
                    del_mssql_table_rows(server, database, samples_table, pk_df=wq_data_sample1[['SiteID', 'Param', 'CollectionTime']])
                to_mssql(wq_data_sample1, server, database, samples_table, dtype={'Value': String})
                print(str(len(wq_data_sample1)) + ' samples updated')

            # Check if data is the same - mtypes
            old_mtype_data = rd_sql(server, database, mtypes_table, where_col={'SiteID': hts_sites.tolist()})
            old_mtype_data.rename(columns={'Value': 'old'}, inplace=True)
            mtype_combo = merge(wq_data_mtype, old_mtype_data, on=mtypes_pkeys, how='outer')
            bad_old_mtype = mtype_combo[mtype_combo.Value.isnull()]
            if not bad_old_mtype.empty:
                del_mssql_table_rows(server, database, mtypes_table, pk_df=bad_old_mtype[['SiteID', 'MeasurementType', 'Param', 'CollectionTime']])
                mtype_combo = mtype_combo[mtype_combo.Value.notnull()]

            wq_data_mtype1 = mtype_combo[mtype_combo.Value != mtype_combo.old].drop('old', axis=1).copy()

            if not wq_data_mtype1.empty:
                if not mtypes_bool1:
                    del_mssql_table_rows(server, database, mtypes_table, pk_df=wq_data_mtype1[['SiteID', 'MeasurementType', 'Param', 'CollectionTime']])
                to_mssql(wq_data_mtype1, server, database, mtypes_table)
                print(str(len(wq_data_mtype1)) + ' measurements updated')

    ### log
    end_run = time()
    mins1 = round((end_run - start_run)/60, 2)
    print(str(mins1) + ' min')
    log1 = DataFrame([[today1, server, 'Python WQ upload', '1', 'INFO', 'python', 'completed successfully']], columns=['Date', 'Server', 'Application', 'Thread', 'Level', 'Logger', 'Message'])
    to_mssql(log1, log_server, log_database, log_table)
    print('complete')

except Exception as err:
    err1 = err
    print(err1)
    log2 = DataFrame([[today1, server, 'Python WQ upload', '1', 'INFO', 'python', str(err1)]], columns=['Date', 'Server', 'Application', 'Thread', 'Level', 'Logger', 'Message'])
    to_mssql(log2, log_server, log_database, log_table)
    print('fail')


###########################################
### Testing

#sites_all[(sites_all.SiteID.str.contains('SQ20104', case=False)) & (sites_all.MeasurementType.str.contains('cond', case=False))].sort_values('MeasurementType')
#


















