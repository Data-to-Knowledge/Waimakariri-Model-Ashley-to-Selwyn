# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/03/2018 9:20 AM
"""

from __future__ import division
from core import env
from nitrate_at_key_receptors import get_stocastic_set, get_n_at_points_nc

if __name__ == '__main__':
    #todo run this shit
    outdir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points")
    nsmc_nums = get_stocastic_set(False)
    get_n_at_points_nc(outdir, nsmc_nums, ucn_var_name='river',
                           ucn_nc_path=r"C:\mh_waimak_model_data\emma_con.nc",
                           cbc_nc_path="C:\mh_waimak_model_data\post_filter1_budget.nc",
                           missing_str_obs='pass')
