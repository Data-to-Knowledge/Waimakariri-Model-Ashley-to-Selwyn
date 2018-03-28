# -*- coding: utf-8 -*-
"""
Created on Mon May 08 09:58:51 2017

@author: MichaelEK

Script to extract time series data from Hydstra and save them to sql tables.

Must be run in a 32bit python!
"""
import os
import pandas as pd
from configparser import ConfigParser
from datetime import date, datetime
from pdsql.mssql import rd_sql, to_mssql, create_mssql_table
from pyhydllp import hyd

#############################################
### Parameters
print('load parameters')

py_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

#py_dir = r'E:\ecan\git\HydroPandas\users\MK\tsdata_exports'

ini1 = ConfigParser()
ini1.read([os.path.join(py_dir, os.path.splitext(__file__)[0] + '.ini')])

link_table = str(ini1.get('Input', 'link_table'))
hydstra_server = str(ini1.get('Input', 'hydstra_server'))
hydstra_database = str(ini1.get('Input', 'hydstra_database'))
ini_path = str(ini1.get('Input', 'ini_path'))
dll_path = str(ini1.get('Input', 'dll_path'))

server = str(ini1.get('Output', 'server'))
database = str(ini1.get('Output', 'database'))
daily_table = str(ini1.get('Output', 'daily_table'))
hourly_table = str(ini1.get('Output', 'hourly_table'))

hydstra_code_sql = {'server': server, 'database': database, 'table': link_table}

today1 = str(date.today())

qual_code_convert = [[10, 11, 18, 20, 21, 30, 31, 40, 50, 60], [600, 600, 520, 500, 500, 400, 400, 300, 200, 100]]
qual_code_dict = dict(zip(qual_code_convert[0], qual_code_convert[1]))
#today1 = date(2017, 12, 22)
#server1 = 'SQL2012DEV01'
#database1 = 'HydstraArchive'
#dtype_dict = {'wtemp': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 1)', 'qual_code': 'INT'}, 'flow': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'precip': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 1)', 'qual_code': 'INT'}, 'swl': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'gwl': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'lakel': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}}

daily_dtype = {'Site': 'VARCHAR(29)', 'FeatureMtypeSourceID': 'int', 'Time': 'date', 'Value': 'float', 'QualityCode': 'int'}
daily_pkeys = ['Site', 'FeatureMtypeSourceID', 'Time']

hourly_dtype = {'Site': 'VARCHAR(29)', 'FeatureMtypeSourceID': 'int', 'Time': 'datetime', 'Value': 'float', 'QualityCode': 'int'}
hourly_pkeys = ['Site', 'FeatureMtypeSourceID', 'Time']

sql_daily = {'server': server, 'database': database, 'table': daily_table}
sql_hourly = {'server': server, 'database': database, 'table': hourly_table}

grp_dict = {'daily': {'interval': 'day', 'export': sql_daily}, 'hourly': {'interval': 'hour', 'export': sql_hourly}}

log_server = str(ini1.get('Output', 'log_server'))
log_database = str(ini1.get('Output', 'log_database'))
log_table = str(ini1.get('Output', 'log_table'))

max_date_stmt = "select max(Time) from " + log_table + " where HydroTable='" + daily_table + "' and RunResult='pass'"

#############################################
### Get the SQL-hydstra codes and last date

hydstra_codes = rd_sql(**hydstra_code_sql).set_index('HydstraCode')['FeatureMtypeSourceID'].to_dict()

last_date1 = str(rd_sql(stmt=max_date_stmt, **sql_daily).loc[0][0].date())
#last_date1 = '2017-12-28'
#today2 = '2017-12-28 00:00:00'

print('Last sucessful date is ' + last_date1)

###############################################
### Create tables if they don't exist

daily_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + daily_table + "', 'U')").loc[0][0] is None

if daily_bool1:
        create_mssql_table(server, database, daily_table, dtype_dict=daily_dtype, primary_keys=daily_pkeys)

hourly_bool1 = rd_sql(server, database, stmt="select OBJECT_ID('" + hourly_table + "', 'U')").loc[0][0] is None

if hourly_bool1:
        create_mssql_table(server, database, hourly_table, dtype_dict=hourly_dtype, primary_keys=hourly_pkeys)

#############################################
### Iterate through hydstra codes and save to SQL

hyd1 = hyd(ini_path, dll_path)

for j in grp_dict:
    print('Interval: ' + j)
    runstart = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    try:
        for i in hydstra_codes:
            print('HydstraCode ' + str(i))
            s1 = hyd1.get_ts_data_bulk(hydstra_server, hydstra_database, int(i), from_mod_date=last_date1, to_mod_date=today1, code_convert=hydstra_codes, qual_code_convert=qual_code_dict, **grp_dict[j])

        runtime = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        log1 = pd.DataFrame([[runstart, grp_dict[j]['export']['table'], 'pass', 'all good', last_date1, runtime]], columns=['Time', 'HydroTable', 'RunResult', 'Comment', 'FromTime', 'RunTimeEnd'])
        to_mssql(log1, log_server, log_database, log_table)

    except Exception as err:
        err1 = err
        print(err1)
        log2 = pd.DataFrame([[runstart, grp_dict[j]['export']['table'], 'fail', str(err1), last_date1]], columns=['Time', 'HydroTable', 'RunResult', 'Comment', 'FromTime'])
        to_mssql(log2, log_server, log_database, log_table)



#t7 = read_hdf(join(base_dir, hydstra_code_dict[i] + '_' + str(today1) + '.h5'))


