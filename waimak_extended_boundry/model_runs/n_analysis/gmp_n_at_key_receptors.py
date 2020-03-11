"""
Author: matth
Date Created: 20/04/2018 6:54 AM
"""

from __future__ import division
import env
import netCDF4 as nc
from waimak_extended_boundry.model_run_tools.n_analysis_support.nitrate_at_key_receptors import get_n_at_points_nc

if __name__ == '__main__':
    outdir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\gmp_n_results_at_points")
    temp = nc.Dataset(r"C:\mh_waimak_model_data\GMP_mednload_ucn.nc")
    nsmc_nums = list(temp.variables['nsmc_num'])
    get_n_at_points_nc(outdir, nsmc_nums, ucn_var_name='mednload',
                           ucn_nc_path=r"C:\mh_waimak_model_data\GMP_mednload_ucn.nc",
                           cbc_nc_path=r"C:\mh_waimak_model_data\GMP_cbc.nc",
                           missing_str_obs='pass')

