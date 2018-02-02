"""
Author: matth
Date Created: 17/05/2017 8:55 AM
"""

from __future__ import division
from core import env
import flopy
import pandas as pd
import numpy as np
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, temp_file_dir

ccc = pd.read_table(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\Water supply wells\SelectedCCCwells.txt")
wdc = pd.read_table(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\Water supply wells\SelectedWDC wells.txt")
all_wells = mt.get_all_well_data()

well_list = None # add chch wells

output_wells_ccc = all_wells[np.in1d(all_wells.index,ccc)]
output_wells_wdc = all_wells[np.in1d(all_wells.index,wdc)]


# import concentration data

ucn = flopy.utils.UcnFile('{}/first_mt3d2/cmp_layer/MT3D001.UCN'.format(base_mod_dir))
con_data = ucn.get_data(ucn.get_kstpkper()[-1])
con_data[np.isclose(con_data,1e30)] = np.nan

for well in output_wells_ccc.index:
    idx = (output_wells_ccc.loc[well,'layer'], output_wells_ccc.loc[well,'row'],output_wells_ccc.loc[well,'col'])
    output_wells_ccc.loc[well,'conc'] = con_data[idx]

output_wells_ccc.to_csv("{}/CHCH_CMP_con_first_look.csv".format(temp_file_dir))

for well in output_wells_wdc.index:
    idx = (output_wells_wdc.loc[well,'layer'], output_wells_wdc.loc[well,'row'],output_wells_wdc.loc[well,'col'])
    output_wells_wdc.loc[well,'conc'] = con_data[idx]

output_wells_wdc.to_csv("{}/WDC_CMP_con_first_look.csv".format(temp_file_dir))