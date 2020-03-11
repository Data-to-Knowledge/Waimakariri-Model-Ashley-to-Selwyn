# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/09/2018 8:53 AM
"""

from __future__ import division
import netCDF4 as nc

def ncdump(ncdataset, outpath):
    lines = []
    # overview
    lines.append('#### netcdf overview ####\n')
    lines.append(str(ncdataset.__str__()))
    lines.append('\n \n \n')


    #variables
    vars = ncdataset.variables.keys()
    vars.sort()
    for var in vars:
        lines.append('#####  {}  #####\n'.format(var))
        lines.append(str(ncdataset.variables[var].__str__()))
        lines.append('\n \n \n')

    #dimension values
    dims = ncdataset.dimensions.keys()
    dims.sort()
    for dim in dims:
        lines.append('#### dimension {} values ####\n'.format(dim))
        lines.append(str(ncdataset.variables[dim][:]))
        lines.append('\n \n \n')

    with open(outpath,'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    data = nc.Dataset(r"K:\mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc")
    ncdump(data,r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\modelling_reports\ashely_waimakarriri_model_build\params_netcdf_overview.txt")