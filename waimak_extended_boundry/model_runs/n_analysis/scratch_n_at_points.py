"""
Author: matth
Date Created: 3/04/2018 3:40 PM
"""

from __future__ import division
from waimak_extended_boundry import \
    calculate_con_from_netcdf_well
from waimak_extended_boundry import \
    get_stocastic_set

if __name__ == '__main__':
    sites = ['L35/0349',
             'L35/0850',
             'L35/1195',
             'M34/0223',
             'M35/0132',
             'M35/0698',
             'M35/11936',
             'M35/4757',
             'M35/4795',
             'M35/5440',
             'M35/5869',
             'M35/6295',
             'M35/7878',
             'M35/8479',
             'M35/8567',
             'M35/8744',
             'M35/4682',
             ]
    nsmc_nums = get_stocastic_set(False)
    ucn_nc_path =r"C:\mh_waimak_model_data\mednload_ucn.nc"
    ucn_var_name = 'mednload'
    all_well_data = calculate_con_from_netcdf_well(nsmc_nums, ucn_nc_path,
                                                   ucn_var_name, sites)

    print('done')
