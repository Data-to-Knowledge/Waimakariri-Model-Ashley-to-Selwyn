# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/01/2018 2:24 PM
"""

from __future__ import division
import os
import netCDF4 as nc
import numpy as np
import pandas as pd
from waimak_extended_boundry.model_run_tools.data_extraction.data_from_streams import \
    _get_sw_samp_pts_dict, get_samp_points_df, _get_flux_flow_arrays
from waimak_extended_boundry.extended_boundry_model_tools import smt
from waimak_extended_boundry.model_run_tools.data_extraction.data_at_wells import get_well_positions
from warnings import warn

def calculate_con_from_netcdf_str(nsmc_nums, ucn_nc_path, ucn_var_name, cbc_nc_path, sites,
                                  outpath=None, missing_str_obs='raise'):
    """
    create a dictionary (keys = runtype) of dataframes(index=nsmc realisation, columns = sites) for each runtype in
    :param nsmc_nums: the nsmc number to use or 'all'
    :param ucn_nc_path: path to netcdf file created by make_ucn_netcdf
    :param ucn_var_name: name of the concentration variable in the ucn netcdf file
    :param cbc_nc_path: path to the cell budget file netcdf (needed for any potential drain sites and joint sfr drn)
    :param sites: list of sites (stream sites)
    :param outpath: optional, if None endmember mixing is not saved, else the directory to save 1 csv of the data
    :param missing_str_obs: action upon missing str_obs, raise Keyerror, warn, or pass
    :return:
    """
    # some checks on sites

    ucn_nc_file = nc.Dataset(ucn_nc_path)
    cbc_nc_file = nc.Dataset(cbc_nc_path)

    filter_nums = np.array(cbc_nc_file.variables['nsmc_num'])
    emma_nums = np.array(ucn_nc_file.variables['nsmc_num'])
    if nsmc_nums == 'all':
        ucn_num_idx = np.ones(emma_nums.shape).astype(bool)
        nsmc_filter_idx = np.in1d(filter_nums, emma_nums)
        nsmc_nums = emma_nums
    else:
        if len(set(nsmc_nums))>len(set(emma_nums)):
            nsmc_nums = emma_nums
        nsmc_nums = np.atleast_1d(nsmc_nums)
        ucn_num_idx = np.in1d(emma_nums, nsmc_nums)
        nsmc_filter_idx = np.in1d(filter_nums, nsmc_nums)

    outdata = pd.DataFrame(index=nsmc_nums, columns=sites)
    sw_samp_dict = _get_sw_samp_pts_dict()
    sw_samp_df = get_samp_points_df()
    for site in sites:
        # set up temp variables
        drn_fraction = None
        drn_flux = None
        sfr_fraction = None
        sfr_flux = None

        # get the surface water feature concentrations
        drn_idx, sfr_idx = _get_flux_flow_arrays(site, sw_samp_dict, sw_samp_df)

        if sfr_idx is not None:
            assert sfr_idx.sum() == 1, 'only single boolean true sfr_idxs used, site:{}'.format(site)
            r, c = smt.model_where(sfr_idx)[0]
            try:
                sfr_fraction = np.array(
                    ucn_nc_file.variables['sobs_{}'.format(ucn_var_name)][ucn_num_idx, r, c])
                sfr_flux = np.array(cbc_nc_file.variables['streamflow out'][nsmc_filter_idx, r, c])
            except KeyError as val:
                if missing_str_obs=='raise':
                    raise KeyError(val)
                elif missing_str_obs =='warn':
                    sfr_fraction = np.nan
                    sfr_flux = np.nan
                    warn(val)
                elif missing_str_obs =='pass':
                    sfr_fraction = np.nan
                    sfr_flux = np.nan
                    pass
        if drn_idx is not None:
            # get drain concentration and flow
            drn_con = np.concatenate([np.array(ucn_nc_file.variables[ucn_var_name][ucn_num_idx, 0, r, c])[:, np.newaxis]
                                      for r, c in smt.model_where(drn_idx)], axis=1)
            drn_flow = np.concatenate(
                [np.array(cbc_nc_file.variables['drains'][nsmc_filter_idx, r, c])[:, np.newaxis]
                 for r, c in smt.model_where(drn_idx)], axis=1)
            assert drn_con.shape == drn_flow.shape, 'drn_con and drn_flow, must be the same size'
            # convert to drain concentration across all drain cells
            drn_fraction = ((drn_con * drn_flow).sum(axis=1) / drn_flow.sum(axis=1))
            drn_flux = drn_flow.sum(axis=1)*-1

        if drn_fraction is not None and sfr_fraction is not None:
            outcon = (sfr_fraction * sfr_flux + drn_fraction * drn_flux) / (drn_flux + sfr_flux)
        elif drn_fraction is not None:
            outcon = drn_fraction
        elif sfr_fraction is not None:
            outcon = sfr_fraction
        else:
            raise ValueError('one of drn_fraction, sfr_fraction should be not None')
        outdata.loc[:, site] = outcon

    outdata.index.name='nsmc_num'
    if outpath is not None:
        # save the files
        if not os.path.exists(os.path.dirname(outpath)):
            os.makedirs(os.path.dirname(outpath))
        outdata.to_csv(outpath)

    return outdata


def calculate_con_from_netcdf_well(nsmc_nums, ucn_nc_path, ucn_var_name, sites, outpath=None):
    """
    create a dictionary (keys = runtype) of dataframes(index=nsmc realisation, columns = sites) for each runtype in
    :param nsmc_nums: the netcdf numbers to calculate values for or all
    :param ucn_nc_path: path to netcdf file created by make_ucn_netcdf
    :param ucn_var_name: name of the concentration variable in the ucn netcdf file
    :param sites: list of sites (well_nums)
    :param outpath: optional, if None endmember mixing is not saved, else the directory to save 1 csv of the data
    :return:
    """
    ucn_nc_file = nc.Dataset(ucn_nc_path)

    emma_nums = np.array(ucn_nc_file.variables['nsmc_num'])
    if nsmc_nums == 'all':
        num_idx = np.ones(emma_nums.shape).astype(bool)
        nsmc_nums = emma_nums
    else:
        if len(set(nsmc_nums))>len(set(emma_nums)):
            nsmc_nums = emma_nums
        nsmc_nums = np.atleast_1d(nsmc_nums)
        num_idx = np.in1d(emma_nums, nsmc_nums)

    well_locs = get_well_positions(sites)

    outdata = pd.DataFrame(index=nsmc_nums, columns=well_locs.index)
    for site in well_locs.index:
        # get the well concentrations
        l, r, c = well_locs.loc[site,['k','i','j']].astype(int)

        well_fraction = np.array(ucn_nc_file.variables[ucn_var_name][num_idx, l, r, c])

        outdata.loc[:, site] = well_fraction

    outdata.index.name='nsmc_num'
    if outpath is not None:
        # save the files
        if not os.path.exists(os.path.dirname(outpath)):
            os.makedirs(os.path.dirname(outpath))
        outdata.to_csv(outpath)

    return outdata


if __name__ == '__main__':
    testtype=1
    if testtype==0:
        test = calculate_con_from_netcdf_str(nsmc_nums='all',
                                             ucn_nc_path=r"T:\Temp\temp_gw_files\mednload_ucn.nc",
                                             ucn_var_name='mednload',
                                             cbc_nc_path=r"K:\mh_modeling\netcdfs_of_key_modeling_data\post_filter1_cell_budgets.nc",
                                             sites=['kaiapoi_harpers_s', 'kaiapoi_heywards'],
                                             outpath=None)
        print(test.describe())
    elif testtype==1:
        test = calculate_con_from_netcdf_well(nsmc_nums = 7,
                                       ucn_nc_path=r"T:\Temp\temp_gw_files\mednload_ucn.nc",
                                       ucn_var_name='mednload',
                                       sites=['M35/1003', 'M35/6295','BS28/5004'])
        print(test)