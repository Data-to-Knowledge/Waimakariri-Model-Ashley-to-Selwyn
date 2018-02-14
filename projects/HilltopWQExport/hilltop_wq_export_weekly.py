# -*- coding: utf-8 -*-
"""
Created on Tue Oct 03 08:51:46 2017

@author: MichaelEK
"""

import sys
from os import path, getcwd

try:

    ######################################
    ### Change pythonpath and add init

    base_py_path = r'D:\Executables\PythonScripts'
    proj_name = 'HilltopWQExport'
    py_path = path.join(base_py_path, proj_name)

    sys.path.append(py_path)

    init1 = path.join(py_path, '__init__.py')
    fh = open(init1,'a+')

    #######################################
    ### Start the work

    from configparser import ConfigParser
    import pandas as pd
    from datetime import datetime
    from core.ecan_io.hilltop import rd_ht_wq_data
    from core.ecan_io.mssql import to_mssql, create_mssql_table, del_mssql_table_rows, rd_sql

    #############################################
    ### Parameters

    py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))

#    py_dir = r'E:\ecan\git\Ecan.Science.Python.Base\projects\HilltopWQExport'

    ini1 = ConfigParser()
    ini1.read([path.join(py_dir, path.splitext(__file__)[0] + '.ini')])

#    ini1.read([path.join(py_dir, 'hilltop_wq_export_weekly' + '.ini')])

    ## Import
    hts_gw = str(ini1.get('Input', 'hts_gw'))
    hts_sw = str(ini1.get('Input', 'hts_sw'))

    sample_params_list = ['Project', 'Cost Code', 'Technician', 'Sample ID', 'Sample Comment', 'Field Comment', 'Sample Appearance', 'Sample Colour', 'Sample Odour', 'Water Colour', 'Water Clarity']
    mtype_params_list = ['Lab Method', 'Lab Name']

    ## Export

    sites_dtype = {'SiteID': 'VARCHAR(19)', 'MeasurementType': 'VARCHAR(59)', 'Units': 'VARCHAR(29)', 'FromDate': 'DATETIME', 'ToDate': 'DATETIME'}

    samples_dtype = {'SiteID': 'VARCHAR(19)', 'Param': 'VARCHAR(59)', 'Value': 'VARCHAR(255)', 'CollectionTime': 'DATETIME'}

    mtypes_dtype = {'SiteID': 'VARCHAR(19)', 'MeasurementType': 'VARCHAR(59)', 'CollectionTime': 'DATETIME', 'Param': 'VARCHAR(59)', 'Value': 'VARCHAR(255)'}

    col_name_convert = {'site': 'SiteID', 'mtype': 'MeasurementType', 'unit': 'Units', 'start_date': 'FromDate', 'end_date': 'ToDate', 'time': 'CollectionTime'}

    server = str(ini1.get('Output', 'server'))
    database = str(ini1.get('Output', 'database'))
    sites_table = str(ini1.get('Output', 'sites_table'))
    mtypes_table = str(ini1.get('Output', 'mtypes_table'))
    samples_table = str(ini1.get('Output', 'samples_table'))

    sites_pkeys = ['SiteID', 'MeasurementType']
    samples_pkeys = ['SiteID', 'Param', 'CollectionTime']
    mtypes_pkeys = ['SiteID', 'MeasurementType', 'Param', 'CollectionTime']

#    sql_log = {'server': 'SQL2014DEV01', 'database': 'Apps_Trace_Log', 'table': 'Log'}
    log_server = str(ini1.get('Output', 'log_server'))
    log_database = str(ini1.get('Output', 'log_database'))
    log_table = str(ini1.get('Output', 'log_table'))

    ##############################################
    ### Determine lastest date and if tables exist

    today1 = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    sites_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + sites_table + "', 'U')").loc[0][0] is None
    samples_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + samples_table + "', 'U')").loc[0][0] is None
    mtypes_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + mtypes_table + "', 'U')").loc[0][0] is None

    ###############################################
    ### Extract data

    mtype_list = list(mtype_params_list)
    mtype_list.append('data')

    ## GW
    gw_wq_data, gw_sites_info = rd_ht_wq_data(hts_gw, output_site_data=True, sample_params=sample_params_list, mtype_params=mtype_params_list)

    gw_sites_info1 = gw_sites_info.drop(['data_source', 'divisor'], axis=1)

    gw_wq_data1 = gw_wq_data.melt(['site', 'mtype', 'time'], var_name='Param', value_name='Value')

    gw_wq_data_mtype = gw_wq_data1[gw_wq_data1['Param'].isin(mtype_list)]

    gw_wq_data_sample = gw_wq_data1[gw_wq_data1['Param'].isin(sample_params_list)]
    gw_wq_data_sample = gw_wq_data_sample.drop('mtype', axis=1).drop_duplicates(['site', 'Param', 'time'])
    gw_wq_data_sample['Value'] = gw_wq_data_sample['Value'].str.strip()
    gw_wq_data_sample = gw_wq_data_sample[~((gw_wq_data_sample['Value'] == '') | (gw_wq_data_sample['Value'] == '..'))]

    # Rename columns
    gw_sites_info1.rename(columns=col_name_convert, inplace=True)
    gw_wq_data_sample.rename(columns=col_name_convert, inplace=True)
    gw_wq_data_mtype.rename(columns=col_name_convert, inplace=True)

    # Remove duplicate primary keys
    gw_sites_info2 = gw_sites_info1.drop_duplicates(sites_pkeys).copy()
    gw_wq_data_sample2 = gw_wq_data_sample.drop_duplicates(samples_pkeys).copy()
    gw_wq_data_mtype2 = gw_wq_data_mtype.drop_duplicates(mtypes_pkeys).copy()

    ## SW
    sw_wq_data, sw_sites_info = rd_ht_wq_data(hts_sw, output_site_data=True, sample_params=sample_params_list, mtype_params=mtype_params_list)
    sw_sites_info1 = sw_sites_info.drop(['data_source', 'divisor'], axis=1)

    sw_wq_data1 = sw_wq_data.melt(['site', 'mtype', 'time'], var_name='Param', value_name='Value')
    sw_wq_data_mtype = sw_wq_data1[sw_wq_data1['Param'].isin(mtype_list)]

    sw_wq_data_sample = sw_wq_data1[sw_wq_data1['Param'].isin(sample_params_list)]
    sw_wq_data_sample = sw_wq_data_sample.drop('mtype', axis=1).drop_duplicates(['site', 'Param', 'time'])
    sw_wq_data_sample['Value'] = sw_wq_data_sample['Value'].str.strip()
    sw_wq_data_sample = sw_wq_data_sample[~((sw_wq_data_sample['Value'] == '') | (sw_wq_data_sample['Value'] == '..'))]

    # Rename columns
    sw_sites_info1.rename(columns=col_name_convert, inplace=True)
    sw_wq_data_sample.rename(columns=col_name_convert, inplace=True)
    sw_wq_data_mtype.rename(columns=col_name_convert, inplace=True)

    # Remove duplicate primary keys
    sw_sites_info2 = sw_sites_info1.drop_duplicates(sites_pkeys).copy()
    sw_wq_data_sample2 = sw_wq_data_sample.drop_duplicates(samples_pkeys).copy()
    sw_wq_data_mtype2 = sw_wq_data_mtype.drop_duplicates(mtypes_pkeys).copy()

    ###############################################
    ### Dump on sql database

    ## Combine datasets

    sites_info2 = pd.concat([gw_sites_info2, sw_sites_info2])
    wq_data_sample2 = pd.concat([gw_wq_data_sample2, sw_wq_data_sample2])
    wq_data_mtype2 = pd.concat([gw_wq_data_mtype2, sw_wq_data_mtype2])

    # Remove duplicate primary keys
    sites_info2 = sites_info2.drop_duplicates(sites_pkeys).dropna(subset=sites_pkeys)
    wq_data_sample2 = wq_data_sample2.drop_duplicates(samples_pkeys).dropna(subset=samples_pkeys)
    wq_data_mtype2 = wq_data_mtype2.drop_duplicates(mtypes_pkeys).dropna(subset=mtypes_pkeys)

    wq_data_sample2['CollectionTime'] = pd.to_datetime(wq_data_sample2['CollectionTime'])
    wq_data_mtype2['CollectionTime'] = pd.to_datetime(wq_data_mtype2['CollectionTime'])

    ## SQL upload

    if sites_bool1:
        create_mssql_table(server, database, sites_table, dtype_dict=sites_dtype, primary_keys=sites_pkeys)
    else:
        del_mssql_table_rows(server, database, sites_table)
    to_mssql(sites_info2, server, database, sites_table)

    if samples_bool1:
        create_mssql_table(server, database, samples_table, dtype_dict=samples_dtype, primary_keys=samples_pkeys)
    else:
        del_mssql_table_rows(server, database, samples_table)
    to_mssql(wq_data_sample2, server, database, samples_table)

    if mtypes_bool1:
        create_mssql_table(server, database, mtypes_table, dtype_dict=mtypes_dtype, primary_keys=mtypes_pkeys)
    else:
        del_mssql_table_rows(server, database, mtypes_table)
    to_mssql(wq_data_mtype2, server, database, mtypes_table)

    ## log
    log1 = pd.DataFrame([[today1, server, 'Python WQ upload', '1', 'INFO', 'python', 'completed successfully']], columns=['Date', 'Server', 'Application', 'Thread', 'Level', 'Logger', 'Message'])
    to_mssql(log1, log_server, log_database, log_table)
    print('complete')

except Exception as err:
    err1 = err
    print(err1)
    log2 = pd.DataFrame([[today1, server, 'Python WQ upload', '1', 'INFO', 'python', str(err1)]], columns=['Date', 'Server', 'Application', 'Thread', 'Level', 'Logger', 'Message'])
    to_mssql(log2, log_server, log_database, log_table)
    print('fail')



