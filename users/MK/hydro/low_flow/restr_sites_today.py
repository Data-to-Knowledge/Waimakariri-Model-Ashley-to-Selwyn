# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 09:34:22 2017

@author: MichaelEK
"""
import sys
new_path = r'E:\ecan\git\Ecan.Science.Python.Base'
sys.path.insert(0, new_path)

from os.path import join
from datetime import date
from core.allo_use.ros import low_flow_restr
#from core.misc.misc import save_df

############################################
### Parameters

base_dir = r'C:\onedrive\OneDrive - Environment Canterbury\share'
output_csv1 = 'restr_site.csv'
output_csv2 = 'restr_site_band.csv'

###########################################
### Run function

today1 = str(date.today())

basic, complete = low_flow_restr(from_date=today1, to_date=today1)

basic.to_csv(join(base_dir, output_csv1), index=False)
complete.to_csv(join(base_dir, output_csv2), index=False)



