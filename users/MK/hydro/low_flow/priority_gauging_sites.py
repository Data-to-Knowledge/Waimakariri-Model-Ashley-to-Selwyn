# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 09:34:22 2017

@author: MichaelEK
"""
import sys
new_path = r'E:\ecan\git\Ecan.Science.Python.Base'
sys.path.insert(0, new_path)

from os.path import join
from core.allo_use.ros import priority_gaugings
#from core.misc.misc import save_df

############################################
### Parameters

base_dir = r'C:\onedrive\OneDrive - Environment Canterbury\share'
output_csv1 = 'lowflow_gauging_sites.csv'

###########################################
### Run function

num_previous_months = 2

priority = priority_gaugings(num_previous_months)

priority.to_csv(join(base_dir, output_csv1))



