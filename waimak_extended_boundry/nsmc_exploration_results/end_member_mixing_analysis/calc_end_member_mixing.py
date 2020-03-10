# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/11/2017 5:00 PM
"""

from __future__ import division
import env
import os
import netCDF4 as nc
import numpy as np
import pandas as pd
from waimak_extended_boundry import \
    get_all_well_row_col
from waimak_extended_boundry import \
    _get_sw_samp_pts_dict, get_samp_points_df, _get_flux_flow_arrays
from waimak_extended_boundry import smt

runtypes = ('coastal', 'inland', 'river')


def calculate_endmember_mixing(ucn_nc_path, cbc_nc_path, sites, outdir=None):
    """
    create a dictionary (keys = runtype) of dataframes(index=nsmc realisation, columns = sites) for each runtype in
    :param ucn_nc_path: path to netcdf file created by make_ucn_netcdf
    :param cbc_nc_path: path to the cell budget file netcdf (needed for any potential drain sites)
    :param sites: dictionary (keys site name values dictionary keys wells, str_drn values
                  list of wells or str_drn names) or pass empty list
                  e.g. {'group1':{'wells':['BW22/0002'], 'str_drn':[]}}
    :param outdir: optional, if None endmember mixing is not saved, else the directory to save 3 csvs of the data
    :return:
    """
    # some checks on sites
    if not isinstance(sites, dict):
        raise ValueError('expected a dictionary for sites see documentation')
    elif not all([isinstance(e, dict) for e in sites.values()]):
        raise ValueError('expected a dictionary of dictionaries, see documentation')
    if not all([sorted(e.keys()) == sorted(['wells', 'str_drn']) for e in sites.values()]):
        raise ValueError('expected only keys of[wells, str_drn] for every entry.  see documentation')

    ucn_nc_file = nc.Dataset(ucn_nc_path)
    cbc_nc_file = nc.Dataset(cbc_nc_path)

    emma_nums = np.array(ucn_nc_file.variables['nsmc_num'])
    filter_nums = np.array(cbc_nc_file.variables['nsmc_num'])
    nsmc_filter_idx = np.in1d(filter_nums, emma_nums)
    if not (emma_nums == filter_nums[nsmc_filter_idx]).all():
        raise ValueError('expected Ucn netcdf to be a sub set and in the same order re nc_nums as cbb')

    missing = set(runtypes) - set(ucn_nc_file.variables.keys())
    if len(missing) > 1:
        raise ValueError('at least two of (coastal, inland, river) '
                         'must be present nc variables: {}'.format(ucn_nc_file.variables.keys()))

    all_wells = get_all_well_row_col()

    outdict = {}
    for runtype in runtypes:
        outdict[runtype] = pd.DataFrame(index=ucn_nc_file.variables['nsmc_num'], columns=sites.keys())
    nsmc_size = ucn_nc_file.dimensions['nsmc_num'].size
    sw_samp_dict = _get_sw_samp_pts_dict()
    sw_samp_df = get_samp_points_df()
    for runtype in (set(runtypes) - missing):
        for site in sites.keys():
            # set up temp variables
            wells = sites[site]['wells']
            str_drns = sites[site]['str_drn']
            well_fraction = []
            drn_fraction = []
            sfr_fraction = []

            if len(wells) != 0:
                # get the well concentrations
                layers, rows, cols = all_wells.loc[wells, ['layer', 'row', 'col']].values.transpose()

                well_fraction = np.concatenate([np.array(ucn_nc_file.variables[runtype][:, l, r, c])[:, np.newaxis]
                                                for l, r, c in zip(layers, rows, cols)], axis=1)
                if well_fraction.shape != (nsmc_size, len(layers)):
                    raise ValueError('weird shape {} expected {}'.format(well_fraction.shape, (nsmc_size, len(layers))))

            if len(str_drns) != 0:
                # get the surface water feature concentrations
                drn_idx = smt.get_empty_model_grid().astype(bool)
                sfr_idx = smt.get_empty_model_grid().astype(bool)
                for str_drn in sites[site]['str_drn']:
                    drn_array, sfr_array = _get_flux_flow_arrays(str_drn, sw_samp_dict, sw_samp_df)
                    if drn_array is not None:
                        drn_idx = drn_idx | drn_array
                    if sfr_array is not None:
                        sfr_idx = sfr_idx | sfr_array

                if sfr_idx.any():
                    raise NotImplementedError('sfr sites not implemented as the data was not captured')

                if drn_idx.any():
                    # get drain concentration and flow
                    drn_con = np.concatenate([np.array(ucn_nc_file.variables[runtype][:, 0, r, c])[:, np.newaxis]
                                              for r, c in smt.model_where(drn_idx)], axis=1)
                    drn_flow = np.concatenate(
                        [np.array(cbc_nc_file.variables['drains'][nsmc_filter_idx, r, c])[:, np.newaxis]
                         for r, c in smt.model_where(drn_idx)], axis=1)
                    # some checks
                    if drn_con.shape != (nsmc_size, len(smt.model_where(drn_idx))):
                        raise ValueError('weird shape for drn con {} expected {}'.format(well_fraction.shape,
                                                                                         (nsmc_size,
                                                                                          len(np.where(drn_idx)[0]))))
                    if drn_flow.shape != (nsmc_size, len(smt.model_where(drn_idx))):
                        raise ValueError('weird shape for drn_flow {} expected {}'.format(well_fraction.shape,
                                                                                          (nsmc_size,
                                                                                           len(np.where(drn_idx)[0]))))
                    # convert to drain concentration across all drain cells
                    drn_fraction = ((drn_con * drn_flow).sum(axis=1) / drn_flow.sum(axis=1))[:, np.newaxis]

            temp_data = np.concatenate([e for e in (well_fraction, drn_fraction, sfr_fraction) if len(e) != 0], axis=1)
            outdict[runtype].loc[:, site] = temp_data.mean(axis=1)

    if len(missing) == 1:
        # calculate the percentage of the missing endmember tracer assuming all sum to 1
        missing_data = outdict[list(missing)[0]]
        missing_data.loc[:] = 1
        for runtype in (set(runtypes) - missing):
            missing_data = missing_data - outdict[runtype]
        outdict[list(missing)[0]] = missing_data

    if outdir is not None:
        # save the files
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for runtype in runtypes:
            outdict[runtype].to_csv(os.path.join(outdir, '{}_end_member_mixing.csv'.format(runtype)))

    return outdict


def end_member_mixing_filter(data_dict, site_limits, method='all_pass', weights=None, phi_threshold=None):
    """
    return a boolean series (index nsmc_num) if passes the end_mixing_filter
    :param data_dict: the output of calculate_endmember_mixing
    :param site_limits: (lower,upper) a tuple dataframes (index sites, columns upper or lower limit of each group)
    :param method: one of: 'all_pass': all sites must be between the specified lower and upper limits
                           'river': the alpine river component must between the specified lower and upper limits
                           'weighted_all_pass': not implemented
                           'weighted_river': not implemented
    :param weights: not used as for weighted
    :param phi_threshold: not used as for weighted
    :return:
    """

    # probably depreciated
    sites = data_dict.values[0].keys()
    nsmc_nums = data_dict.values[0].index
    # input data checks
    if method not in ['all_pass', 'river']:
        raise NotImplementedError('method {} not implemented'.format(method))

    if len(site_limits) != 2:
        raise ValueError('must pass both the upper and lower Dataframes or Series')
    for df in site_limits.values:
        if not isinstance(df, pd.DataFrame):
            raise ValueError('upper and lower must be a Dataframe or Series')
        if set(df.index) != set(sites):
            raise ValueError('upper and lower must be defined for all sites: {}'.format(sites))
        if set(df.keys()) != set(runtypes):
            raise ValueError('upper and lower must be defined for all runtypes: {}'.format(runtypes))

    if method == 'all_pass':
        inlimits = np.zeros((len(nsmc_nums), 3)).astype(bool)
        for i, rt in enumerate(runtypes):
            temp = (data_dict[rt] >= site_limits[0].loc[:, rt]) & (
                data_dict[rt] <= site_limits[0].loc[:, rt])
            inlimits[:, i] = temp
        passed = inlimits.all(axis=1)

    elif method == 'river':
        rt = 'river'
        passed = (data_dict[rt] >= site_limits[0].loc[:, rt]) & (
            data_dict[rt] <= site_limits[0].loc[:, rt])

    else:
        raise ValueError('should not get here')

    passed = pd.Series(index=nsmc_nums, data=passed)
    return passed


if __name__ == '__main__':
    ucn_path = env.gw_met_data("mh_modeling/netcdfs_of_key_modeling_data/emma_unc.nc")
    cbc_path = env.gw_met_data(
        "mh_modeling/netcdfs_of_key_modeling_data/post_filter1_cell_budgets.nc")
    # set up well sites
    site_info = pd.read_excel(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Sites_GroupList.xlsx"))
    site_info.loc[:, 'well_group'] = [e.lower() for e in site_info.Group]
    sites = {}
    for g in set(site_info.well_group):
        if g == 'none':
            continue
        sites[g] = {'wells': list(site_info.loc[site_info.well_group == g, 'Well list']), 'str_drn': []}

    # set up sw sites
    sw_site_names = ['ohoka_island', 'kaiapoi_island', 'courtenay_swaz', 'ohoka_offbutch', 'kaiapoi_neeves',
                     'ohoka_jeffs', 'kaiapoi_heywards', 'kaiapoi_dixons']
    for sw in sw_site_names:
        sites[sw] = {'wells': [], 'str_drn': [sw]}
    base_dir = env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\nsmc_results\emma_filter")
    data_dict = calculate_endmember_mixing(ucn_nc_path=ucn_path, cbc_nc_path=cbc_path, sites=sites,
                                           outdir=base_dir)

