# -*- coding: utf-8 -*-
"""
Created on Wed May 10 10:45:39 2017

@author: MichaelEK
"""


############################################
#### NIWA Climate data

from xarray import open_dataset, open_mfdataset
from os.path import join
from core.misc import rd_dir

base_dir = r'C:\ecan\niwa\RCPpast\BCC-CSM1.1'
nc1 = 'TotalPrecipCorr_VCSN_BCC-CSM1.1_RCPpast_1971_2005_south-island_p05_daily_ECan.nc'

ncs = rd_dir(base_dir, '.nc')

x1 = open_dataset(join(base_dir, ncs[2]))
x1
x1.close()

x0 = open_mfdataset(join(base_dir, '*.nc'))
x0
x0.close()

top_nc1 = r'C:\ecan\niwa\Test_Outputs\streamq_daily_average_1972010200_1985123100_utc_topnet_13036997_strahler3-B1.nc'

x2 = open_dataset(top_nc1)
x2
x2.close()

rcp1 = r'C:\ecan\niwa\RCP2.6\BCC-CSM1.1\MaxTempCorr_VCSN_BCC-CSM1.1_RCP2.6_2006_2120_south-island_p05_daily_ECan.nc'

x3 = open_dataset(rcp1)
x3
x3.close()

rcp2 = r'C:\ecan\niwa\RCP2.6\BCC-CSM1.1\PE_VCSN_BCC-CSM1.1_RCP2.6_2006_2100_south-island_p05_daily_ECan.nc'

x4 = open_dataset(rcp2)
x4
x4.close()

rcp2 = r'C:\ecan\ftp\niwa\MaxTempCorr_VCSN_BCC-CSM1.1_RCP2.6_2006_2120_south-island_p05_daily_ECan.nc'

x4 = open_dataset(rcp2)
x4.close()


######################################################
#### Testing nc files

### Parameters
base_dir = r'C:\ecan\ftp\niwa'
dir1 = ['RCPpast', 'RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
dir2 = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
dir3 = ['GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']

### Reading nc files

dir0 = dir1[2]
for i in dir2:
    path1 = join(base_dir, dir0, i)
    files = rd_dir(path1, 'nc')
    for j in files:
        path2 = join(path1, j)
        x0 = open_dataset(path2)
        print(x0)
        x0.close()


#########################################################
#### Extract zip files in all subdirectories

import zipfile, fnmatch, os

rootPath = r"C:\ecan\ftp\niwa\RCP6.0"
pattern = '*.zip'

for root, dirs, files in os.walk(rootPath):
    for filename in fnmatch.filter(files, pattern):
        print(os.path.join(root, filename))
        zipfile.ZipFile(os.path.join(root, filename)).extractall(root)

##  and test

for i in dir2:
    path1 = join(rootPath, i)
    files = rd_dir(path1, 'nc')
    for j in files:
        path2 = join(path1, j)
        x0 = open_dataset(path2)
        print(x0)
        x0.close()









rootPath = r"C:\ecan\ftp\niwa\RCP8.5"
pattern = '*.zip'

for root, dirs, files in os.walk(rootPath):
    for filename in fnmatch.filter(files, pattern):
        print(os.path.join(root, filename))
        zipfile.ZipFile(os.path.join(root, filename)).extractall(root)

##  and test

for i in dir2:
    path1 = join(rootPath, i)
    files = rd_dir(path1, 'nc')
    for j in files:
        path2 = join(path1, j)
        x0 = open_dataset(path2)
        print(x0)
        x0.close()


### Problem files

f1 = 'C:\ecan\ftp\niwa\RCP6.0\CESM1-CAM5\MSLP_VCSN_CESM1-CAM5_RCP6.0_2006_2120_south-island_p05_daily_ECan.zip'

##############################################
### TopNet

from core.misc import unarchive_dir

### Parameters
base_dir = r'C:\ecan\ftp\niwa\topnet'
dir1 = ['RCPpast', 'RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
dir2 = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
dir3 = ['GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']

### Reading nc files

dir0 = dir1[4]
for i in dir2:
    path1 = join(base_dir, dir0, i)
    files = rd_dir(path1, 'nc')
    for j in files:
        path2 = join(path1, j)
        x0 = open_dataset(path2)
        print(x0)
        x0.close()



path1 = r'I:\niwa_data\topnet\Climate'

unarchive_dir(path1, 'gz', True)

folder = path1
ext = 'gz'

import fnmatch, os

for root, dirs, files in os.walk(folder):
    for filename in fnmatch.filter(files, '*.' + ext):
        print(os.path.join(root, filename))
        os.remove(os.path.join(root, filename))


##########################################
#### New Topnet set

from xarray import open_dataset, open_mfdataset, Dataset, concat
from os import path, walk
from core.misc import rd_dir
import fnmatch

### Parameters
base_dir = r'I:\niwa_data\topnet\Climate'
base_dir = r'I:\niwa_data\topnet\Climate\RCP8.5\HadGEM2-ES'

### Reading nc files

ext = 'nc'

for root, dirs, files in walk(base_dir):
    for filename in fnmatch.filter(files, '*.' + ext):
        if '_snw' not in filename:
            path2 = path.join(root, filename)
            x0 = open_dataset(path2)
#            print(x0)
            if 'x1' in locals():
                x1 = x1.merge(x0)
            else:
                x1 = x0.copy()
#            x0.close()
    if 'x1' in locals():
        print(x1)
        x1.close()
        del x1



j = files[11]









































