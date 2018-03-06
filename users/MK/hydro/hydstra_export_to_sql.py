# -*- coding: utf-8 -*-
"""
Created on Mon May 08 09:58:51 2017

@author: MichaelEK

Script to extract time series data from Hydstra and save them to sql tables.

Must be run in a 32bit python!
"""

from core.ecan_io import rd_sql, write_sql
from core.ecan_io.hydllp import rd_hydstra
from pandas import concat, read_hdf, read_csv, merge, Timestamp, to_datetime, DataFrame, to_numeric, DateOffset
from os.path import join
from datetime import date
from core.misc.misc import save_df
from time import time

#############################################
### Parameters
base_dir = r'\\fs02\ManagedShares2\Data\Surface Water\shared\base_data'

today1 = date.today()
today1 = date(2017, 12, 22)
#server1 = 'SQL2012DEV01'
#database1 = 'HydstraArchive'
#dtype_dict = {'wtemp': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 1)', 'qual_code': 'INT'}, 'flow': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'precip': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 1)', 'qual_code': 'INT'}, 'swl': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'gwl': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}, 'lakel': {'site': 'VARCHAR(19)', 'time': 'DATE', 'data': 'NUMERIC(10, 3)', 'qual_code': 'INT'}}

interval = 'hour'
export = {'server': 'SQL2012DEV01', 'database': 'Hydro', 'table': 'HydstraTSDataHourly'}
#export = {'server': 'SQL2012DEV01', 'database': 'Hydro', 'table': 'HydstraTSDataDaily'}
hydstra_code_list = [10]

#############################################
### Iterate through hydstra codes and save as hdf files

start1 = time()
for i in hydstra_code_list:
    print('Hydstra code: ' + str(i))
#    s1 = rd_hydstra(i, interval=interval, export=join(base_dir, hydstra_code_dict[i] + '_' + str(today1) + '.h5'), concat_data=False)
    s1 = rd_hydstra(i, interval=interval, export=export, concat_data=False)
#    save_df(s1, join(base_dir, hydstra_code_dict[i] + '_' + str(today1) + '.h5'))

end1 = time()
end1 - start1


#t7 = read_hdf(join(base_dir, hydstra_code_dict[i] + '_' + str(today1) + '.h5'))


