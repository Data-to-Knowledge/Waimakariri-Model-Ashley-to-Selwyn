"""
 Author: Matt Hanson
 Created: 3/3/2020 11:43 AM
 """

# just a test # second test

import netCDF4 as nc
import numpy as np
import glob

if __name__ == '__main__':
    paths = glob.glob(r"D:\Waimakariri_model_input_data\recommended\*.nc")
    for p in paths:
        print('\n\n')
        print(p)
        temp = nc.Dataset(p)
        for k in temp.variables:
            print(k)