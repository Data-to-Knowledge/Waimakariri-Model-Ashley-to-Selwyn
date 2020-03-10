# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 13/11/2017 12:12 PM
"""

from __future__ import division
import env
import numpy as np
import netCDF4 as nc
import flopy_mh as flopy
from waimak_extended_boundry import \
    hds_no_data

# quick check for the netcdf
def spot_check_hds(n):
    data = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\post_filter1_hds.nc"))
    nsmc_nums = data.variables['nsmc_num']
    check_nums = np.random.random_integers(1, len(nsmc_nums) - 2, n)

    for i in check_nums:
        test_num = nsmc_nums[i]
        hds = flopy.utils.HeadFile(r"K:\mh_modeling\data_from_gns\hdsrepo\mf_aw_ex_{}.hds".format(test_num)).get_data(
            (0, 0))
        hds[np.isclose(hds, hds_no_data)] = np.nan

        nc_hds = np.array(data['heads'][i])
        out = (np.isclose(hds, nc_hds) | (np.isnan(hds) & np.isnan(nc_hds))).all()
        print('{}: {}'.format(test_num, out))


def spot_check_cbc_non_flow(n):
    data = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\post_filter1_cell_budgets.nc"))
    nsmc_nums = data.variables['nsmc_num']
    check_nums = np.random.random_integers(1, len(nsmc_nums) - 2, n)

    for var in ['flow right face', 'constant head', 'flow front face', 'flow lower face', 'wells', 'drains',
                'recharge', 'stream leakage']:
        for i in check_nums:
            test_num = nsmc_nums[i]
            cbc = flopy.utils.CellBudgetFile(r"K:\mh_modeling\data_from_gns\cbcrepo\mf_aw_ex_{}.cbc".format(test_num)
                                             ).get_data(kstpkper=(0, 0), text=var, full3D=True)[0]
            if isinstance(cbc, np.ma.MaskedArray):
                cbc = cbc.filled(np.nan)
            cbc = cbc[0] # just first layer

            nc_cbc = np.array(data[var][i])
            if nc_cbc.ndim == 3:
                nc_cbc=nc_cbc[0]
            out = (np.isclose(cbc, nc_cbc) | (np.isnan(cbc) & np.isnan(nc_cbc))).all()
            print('{} {}: {}'.format(var, test_num, out))

def spot_check_cbc_flow(n):
    data = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\post_filter1_cell_budgets.nc"))
    nsmc_nums = data.variables['nsmc_num']
    check_nums = np.random.random_integers(1, len(nsmc_nums) - 2, n)

    for var in ['streamflow out']:
        for i in check_nums:
            test_num = nsmc_nums[i]
            cbc = flopy.utils.CellBudgetFile(r"K:\mh_modeling\data_from_gns\sforepo\mf_aw_ex_{}.sfo".format(test_num)
                                             ).get_data(kstpkper=(0, 0), text=var, full3D=True)[0]
            if isinstance(cbc, np.ma.MaskedArray):
                cbc = cbc.filled(np.nan)
            cbc = cbc[0] # just first layer

            nc_cbc = np.array(data[var][i])
            if nc_cbc.ndim == 3:
                nc_cbc=nc_cbc[0]
            out = (np.isclose(cbc, nc_cbc) | (np.isnan(cbc) & np.isnan(nc_cbc))).all()
            print('{} {}: {}'.format(var, test_num, out))


if __name__ == '__main__':
    spot_check_cbc_flow(5)
