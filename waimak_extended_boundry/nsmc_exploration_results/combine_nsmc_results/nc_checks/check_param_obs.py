# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 10:41 AM
"""

from __future__ import division
import netCDF4 as nc


# quick check for the netcdf
def check_attributes(nc_path, outpath):
    # write out every variable and attribute of the nc file to a text file so I can check it
    f = open(outpath,'w')
    data = nc.Dataset(nc_path)
    variables = data.variables.keys()

    # write general attributes
    gen_atts = data.ncattrs()
    f.write('######################################\n')
    f.write('GENERAL ATTRIBUTES\n')
    f.write('######################################\n\n')

    for at in gen_atts:
        f.write('######################################\n\n')
        f.write('{}\n\n'.format(at))
        f.write('{}\n\n'.format(data.getncattr(at)))

    # write attributes of each variable
    f.write('######################################\n')
    f.write('variable ATTRIBUTES\n')
    f.write('######################################\n\n')
    for var in variables:
        f.write('######################################\n\n')
        f.write('{}\n\n'.format(var))

        atts = data.variables[var].ncattrs()
        for at in atts:
            f.write('{}: {}\n'.format(at,data.variables[var].getncattr(at)))
    f.close()

def check_random_entries():
    # write out data frames with n random entries for each variable
    raise NotImplementedError

if __name__ == '__main__':
    check_attributes(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc",
                     r"C:\Users\MattH\Downloads\param_oobs_att_check.txt")