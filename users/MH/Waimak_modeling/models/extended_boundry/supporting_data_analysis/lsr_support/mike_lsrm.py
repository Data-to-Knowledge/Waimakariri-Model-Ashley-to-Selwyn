# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/10/2017 2:04 PM
"""

from __future__ import division
from core import env
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:42:42 2017

@author: MichaelEK

Example script for running the LSRM.
"""

from core.ecan_io import rd_niwa_data_lsrm
from time import time
from core.ts.gw.lsrm import poly_import, input_processing, lsrm
from core.misc import save_df, rd_dir, get_subdir
from itertools import product
from os import path
from pandas import read_hdf
from collections import OrderedDict
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
start1 = time()

#############################################
#### Parameters

### Reading data
#irr_type_dict = {'server': 'SQL2012PROD05', 'database': 'GIS',
                # 'table': 'AQUALINC_NZTM_IRRIGATED_AREA_20160629', 'column': 'type'} # not useing the 2016 irrigation
irr_type_dict = {'shp': "{}\m_ex_bd_inputs\shp\wai_irr_area_intersect.shp".format(smt.sdp), 'column': 'type'}  # note I have set all irrigation types to Pivot
paw_dict = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'LAND_NZTM_NEWZEALANDFUNDAMENTALSOILS', 'column': 'PAW_MID'}

bound_shp = r'P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\model_grid_domain.shp'


buffer_dis = 10000
from_date = None
to_date = None
paw_ratio = 0.67

include_irr = True

### input data processing
rain_name = 'rain'
pet_name = 'pe'

time_agg = 'W' # Use 'D' for day, 'W' for week, or 'M' for month
agg_ts_fun = 'sum'
grid_res = 1000
crs = 4326
interp_fun = 'linear'
precip_correction = 1.1
min_irr_area_ratio = 0.01

irr_mons = [10, 11, 12, 1, 2, 3, 4]

irr_eff_dict = {'Drip/micro': 1, 'Unknown': 0.8, 'Gun': 0.8, 'Pivot': 0.8, 'K-line/Long lateral': 0.8, 'Rotorainer': 0.8, 'Solid set': 0.8, 'Borderdyke': 0.5, 'Linear boom': 0.8, 'Unknown': 0.8, 'Lateral': 0.8, 'Wild flooding': 0.5, 'Side Roll': 0.8}
irr_trig_dict = {'Drip/micro': 0.5, 'Unknown': 0.5, 'Gun': 0.5, 'Pivot': 0.5, 'K-line/Long lateral': 0.5, 'Rotorainer': 0.5, 'Solid set': 0.5, 'Borderdyke': 0.5, 'Linear boom': 0.5, 'Unknown': 0.5, 'Lateral': 0.5, 'Wild flooding': 0.5, 'Side Roll': 0.5}

### Model parameters
A = 6

### Output parameters
output_shp = r'D:\lsrm_results_v2\test\output_test2.shp'
output_dir = r'D:\lsrm_results_v2'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(os.path.dirname(output_shp)):
    os.makedirs(os.path.dirname(output_shp))

with open(os.path.join(output_dir,'READ_ME.txt'),'w') as f:
    f.write('using irrigation data {}'.format(irr_type_dict))


### Iteration parameters
rcppast = r'D:\niwa_data\climate_projections\RCPpast'
rcp4 = r'D:\niwa_data\climate_projections\RCP4.5'
rcp8 = r'D:\niwa_data\climate_projections\RCP8.5'
vcsn = r'\\fileservices02\ManagedShares\Data\VirtualClimate\vcsn_precip_et_2016-06-06.nc'

all_dir = []
rcppast_dir = get_subdir(rcppast, True)
rcp4_dir = get_subdir(rcp4, True)
rcp8_dir = get_subdir(rcp8, True)
all_dir.extend(rcp4_dir)
all_dir.extend(rcp8_dir)
all_dir.extend(rcppast_dir)

base_param_dict = {'no_irr': {'irr_eff': 1, 'include_irr': False}, '80perc': {'irr_eff': 0.8, 'include_irr': True}, '100perc': {'irr_eff': 1, 'include_irr': True}}
#vcsn_param_dict = {'50perc': {'irr_eff': 0.5, 'include_irr': True}, '65perc': {'irr_eff': 0.65, 'include_irr': True}, '80perc': {'irr_eff': 0.8, 'include_irr': True}, '100perc': {'irr_eff': 1, 'include_irr': True}, 'no_irr': {'irr_eff': 1, 'include_irr': False}}
vcsn_param_dict = {'80perc': {'irr_eff': 0.8, 'include_irr': True}, '100perc': {'irr_eff': 1, 'include_irr': True}, 'no_irr': {'irr_eff': 1, 'include_irr': False}}

param_dict = OrderedDict()
#for f in all_dir: # this is for projectsion
#    model = path.split(f)[1]
#    rcp = path.split(path.split(f)[0])[1]
#    nc_dir = f
#    for d in base_param_dict:
#        name = '_'.join([rcp, model, d])
#        dict1 = base_param_dict[d].copy()
#        dict1.update({'nc_dir': nc_dir, 'from_date': None, 'to_date': None})
#        param_dict.update({name: dict1})

for v in vcsn_param_dict: # this is for VCSN
    name = '_'.join(['vcsn', v])
    dict1 = vcsn_param_dict[v].copy()
    dict1.update({'nc_dir': vcsn, 'from_date': '2008-07-01', 'to_date': '2015-06-30'})
    param_dict.update({name: dict1})

irr1, paw1 = poly_import(irr_type_dict, paw_dict, paw_ratio)

##########################################
### Run through many iterations
for it in param_dict:
    nc_dir = param_dict[it]['nc_dir']
    include_irr = param_dict[it]['include_irr']
    output_results = path.join(output_dir, it + '.h5')
    irr_eff = param_dict[it]['irr_eff']
    irr_eff_dict = {i: irr_eff for i in irr_eff_dict}
    from_date = param_dict[it]['from_date']
    to_date = param_dict[it]['to_date']

    print('Start run for ' + it)

    ###########################################
    ### Extract data

    print('Read in the input data')

    precip_et = rd_niwa_data_lsrm(bound_shp, nc_dir, buffer_dis=buffer_dis, from_date=from_date, to_date=to_date)
    ##########################################
    ### Process data
    ## Resample met data data

    print('Process the input data for the LSR model')

    model_var, sites_poly = input_processing(precip_et, crs, irr1, paw1, bound_shp, rain_name, pet_name, grid_res, buffer_dis, interp_fun, agg_ts_fun, time_agg, irr_eff_dict, irr_trig_dict, min_irr_area_ratio, irr_mons, precip_correction)

    #########################################################
    #### Run the model

    print('Run the LSR model')

    output1 = lsrm(model_var, A, include_irr=include_irr)

    #######################################################
    #### Output data

    print('Saving data')

    save_df(output1, output_results, index=False)

    del(precip_et, model_var, output1)

sites_poly.reset_index().to_file(output_shp) #this outputs the polygons, so only needs running once don't delete

end1 = time()

diff = end1 - start1




