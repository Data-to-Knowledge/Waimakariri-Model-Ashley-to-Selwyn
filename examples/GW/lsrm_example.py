# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 10:42:42 2017

@author: MichaelEK

Example script for running the LSRM.
"""

from core.ecan_io import rd_niwa_vcsn
from time import time
from core.ts.gw.lsrm import poly_import, input_processing, lsrm
from core.misc import save_df

#############################################
#### Parameters

### Reading data
irr_type_dict = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'AQUALINC_NZTM_IRRIGATED_AREA_20160629', 'column': 'type'}
paw_dict = {'server': 'SQL2012PROD05', 'database': 'GIS', 'table': 'LAND_NZTM_NEWZEALANDFUNDAMENTALSOILS', 'column': 'PAW_MID'}

bound_shp = r'\\fileservices02\ManagedShares\Data\VirtualClimate\examples\waipara.shp'

buffer_dis = 10000
from_date = '1972-01-01'
to_date = '2016-06-06'
paw_ratio = 0.67

### input data processing
rain_name = 'rain'
pet_name = 'pe'

time_agg = 'W' # Use 'D' for day, 'W' for week, or 'M' for month
agg_ts_fun = 'sum'
grid_res = 1000
crs = 4326
interp_fun = 'cubic'
precip_correction = 1.1
min_irr_area_ratio = 0.01

irr_mons = [10, 11, 12, 1, 2, 3, 4]

irr_eff_dict = {'Drip/micro': 1, 'Unknown': 0.8, 'Gun': 0.8, 'Pivot': 0.8, 'K-line/Long lateral': 0.8, 'Rotorainer': 0.8, 'Solid set': 0.8, 'Borderdyke': 0.5, 'Linear boom': 0.8, 'Unknown': 0.8, 'Lateral': 0.8, 'Wild flooding': 0.5, 'Side Roll': 0.8}
irr_trig_dict = {'Drip/micro': 0.7, 'Unknown': 0.5, 'Gun': 0.5, 'Pivot': 0.5, 'K-line/Long lateral': 0.5, 'Rotorainer': 0.5, 'Solid set': 0.5, 'Borderdyke': 0.5, 'Linear boom': 0.5, 'Unknown': 0.5, 'Lateral': 0.5, 'Wild flooding': 0.5, 'Side Roll': 0.5}

### Model parameters
A = 6

### Output parameters
output_results = r'\\fileservices02\ManagedShares\Data\VirtualClimate\examples\output_waipara1.h5'
output_shp = r'\\fileservices02\ManagedShares\Data\VirtualClimate\examples\output_waipara1.shp'

###########################################
### Extract data

print('Read in the input data')

irr1, paw1 = poly_import(irr_type_dict, paw_dict, paw_ratio)

precip_et = rd_niwa_vcsn(mtypes=['precip', 'PET'], sites=bound_shp, buffer_dis=buffer_dis, from_date=from_date, to_date=to_date)

##########################################
### Process data
## Resample met data data

print('Process the input data for the LSR model')

model_var, sites_poly = input_processing(precip_et, crs, irr1, paw1, bound_shp, rain_name, pet_name, grid_res, buffer_dis, interp_fun, agg_ts_fun, time_agg, irr_eff_dict, irr_trig_dict, min_irr_area_ratio, irr_mons, precip_correction)

#########################################################
#### Run the model

print('Run the LSR model')

output1 = lsrm(model_var, A)

#######################################################
#### Output data

print('Saving data')

save_df(output1, output_results, index=False)
sites_poly.reset_index().to_file(output_shp)






