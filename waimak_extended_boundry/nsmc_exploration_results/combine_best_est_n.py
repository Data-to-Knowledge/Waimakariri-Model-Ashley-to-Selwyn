# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/11/2017 10:49 AM
"""

from __future__ import division
import env
from waimak_extended_boundry import make_ucn_netcd, emma_nsmc_numbers
import os
from future.builtins import input

def netcdf_best_est_n(nc_path,zlib):
    """
    make the netcdf for the End member mixing analysis
    :param nc_path: path to save the netcdf
    :return:
    """
    vars = ['mednload']
    vars_base_paths = [
        env.gw_met_data("mh_modeling/data_from_gns/median_n/MedNload_ucnrepo"),
    ]
    stobs_base_paths = [env.gw_met_data("mh_modeling\data_from_gns\median_n\MedNload_sobsrepo")]
    stobs_paths = {}
    ucn_paths = {}
    for var, bp, sfobp in zip(vars, vars_base_paths, stobs_base_paths):
        temp_paths = []
        temp_sfo_paths = []
        for i in emma_nsmc_numbers:
            if i > 0:
                temp_paths.append(os.path.join(bp, 'mt_aw_ex_{}_{}.ucn'.format(var, i)))
                temp_sfo_paths.append(os.path.join(sfobp, 'mt_aw_ex_{}_{}.sobs'.format(var, i)))
            elif i == -1:
                temp_paths.append(os.path.join(bp, 'mt_aw_ex_{}_philow.ucn'.format(var)))
                temp_sfo_paths.append(os.path.join(sfobp, 'mt_aw_ex_{}_philow.sobs'.format(var)))
            elif i == -2:
                temp_paths.append(os.path.join(bp, 'mt_aw_ex_{}_phiupper.ucn'.format(var)))
                temp_sfo_paths.append(os.path.join(sfobp, 'mt_aw_ex_{}_phiupper.sobs'.format(var)))
            else:
                raise ValueError('unexpected number in emma_nsmc_nums')
        ucn_paths[var] = temp_paths
        stobs_paths[var] = temp_sfo_paths

    description = """ best estimate nitrate loads at steady state"""

    make_ucn_netcd(nsmc_nums=emma_nsmc_numbers, ucn_paths=ucn_paths, units='g/m3',
                   description=description, nc_path=nc_path, sobs=stobs_paths,zlib=zlib)


if __name__ == '__main__':
    cont = input('are you sure you want to re-run make mednload UcN netcdfs it will overwrite and takes some time y/n')
    if cont != 'y':
        raise ValueError('user interuppted process to prevent overwrite')

    #netcdf_best_est_n(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\mednload_unc.nc"),zlib=True)
    netcdf_best_est_n(r"C:\mh_waimak_model_data\mednload_ucn.nc",zlib=False) # run on gw02
