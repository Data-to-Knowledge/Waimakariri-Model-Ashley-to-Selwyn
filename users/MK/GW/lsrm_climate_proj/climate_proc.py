# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:57:44 2017

@author: Michaelek
"""

from xarray import open_dataset, concat, Dataset
from os import path, walk, makedirs
from pandas import to_datetime
from geopandas import read_file
from core.ecan_io.met import rd_niwa_data_lsrm

###############################################
### Parameters

## Input data

proj_dir = r'D:\niwa_data\climate_projections'

bound_shp = r'P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\model_grid_domain.shp'
bound_shp = r'\\fileservices02\ManagedShares\Data\VirtualClimate\examples\waipara.shp'

bound = read_file(bound_shp)

ds6 = rd_niwa_data_lsrm(bound_shp)





































###################################################
### Testing

nc1 = r'D:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1\TotalPrecipCorr_VCSN_BCC-CSM1.1_RCP2.6_2006_2120_south-island_p05_daily_ECan.nc'
nc2 = r'D:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1\PE_VCSN_BCC-CSM1.1_RCP2.6_2006_2100_south-island_p05_daily_ECan.nc'

ds1 = open_dataset(nc1, drop_variables=['elevation', 'time_bounds', 'time_bnds'])
ds2 = open_dataset(nc2, drop_variables=['elevation', 'time_bounds', 'time_bnds'])

nc_dir = r'D:\niwa_data\climate_projections\RCP2.6\BCC-CSM1.1'






























