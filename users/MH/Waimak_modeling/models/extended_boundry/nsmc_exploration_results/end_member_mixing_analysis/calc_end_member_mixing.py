# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/11/2017 5:00 PM
"""

from __future__ import division
from core import env
import os
import netCDF4 as nc
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_from_streams import \
    _get_sw_samp_pts_dict, get_samp_points_df, _get_flux_flow_arrays
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

runtypes = ('coastal', 'inland', 'river')


def calculate_endmember_mixing(nc_path, sites, outdir=None):
    """
    create a dictionary (keys = runtype) of dataframes(index=nsmc realisation, columns = sites) for each runtype in
    :param nc_path: path to netcdf file created by make_ucn_netcdf
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

    nc_file = nc.Dataset(nc_path)

    missing = set(runtypes) - set(nc_file.variables.keys())
    if len(missing) > 1:
        raise ValueError('at least two of (coastal, inland, river) '
                         'must be present nc variables: {}'.format(nc_file.variables.keys()))

    all_wells = get_all_well_row_col()

    outdict = {}
    for runtype in runtypes:
        outdict[runtype] = pd.DataFrame(index=nc_file.variables['nsmc_num'], columns=sites.keys())
    nsmc_size = nc_file.dimensions['nsmc_num'].size
    sw_samp_dict = _get_sw_samp_pts_dict()
    sw_samp_df = get_samp_points_df()
    for runtype in (set(runtypes) - missing):
        for site in sites.keys:
            wells = sites[site]['wells']
            well_fraction = []
            str_drn_fraction = []
            if len(wells) != 0:
                layers, rows, cols = all_wells.loc[wells, ['layer', 'row', 'col']]  # todo check broadcasting
                well_fraction = np.concatenate([nc_file.variables[runtype][:, l, r, c][:, np.newaxis]
                                                for l, r, c in zip(layers, rows, cols)], axis=1)
                if well_fraction.shape != (nsmc_size, len(layers)):
                    raise ValueError('weird shape {} expected {}'.format(well_fraction.shape, (nsmc_size, len(layers))))

            str_drns = sites[site]['str_drn']
            if len(str_drns) != 0:
                drn_idx = smt.get_empty_model_grid().astype(bool)
                for str_drn in sites[site]['str_drn']:
                    drn_array, sfr_array = _get_flux_flow_arrays(str_drn, sw_samp_dict, sw_samp_df)

                    if sfr_array.any():
                        raise NotImplementedError('sfr sites not implemented as they are hard...')
                    drn_idx = drn_idx | drn_array
                str_drn_fraction = np.concatenate([nc_file.variables[runtype][:, 0, r, c][:, np.newaxis]
                                                   for r, c in zip(np.where(drn_idx))], axis=1)
                if str_drn_fraction.shape != (nsmc_size, len(np.where(drn_idx)[0])):
                    raise ValueError('weird shape {} expected {}'.format(well_fraction.shape,
                                                                         (nsmc_size, len(np.where(drn_idx)[0]))))

            temp_data = np.concatenate((well_fraction, str_drn_fraction), axis=1)
            outdict[runtype].loc[:, site] = temp_data.mean(axis=1)

    if len(missing) == 1:
        # calculate the percentage of the missing well assuming all sum to 1
        missing_data = outdict[list(missing)[0]]
        missing_data.loc[:] = 1
        for runtype in (set(runtypes) - missing):
            missing_data = missing_data - outdict[runtype]

    if outdir is not None:
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
            data_dict[rt] <= site_limits[0].loc[:, rt])  # todo check
            inlimits[:, i] = temp
        passed = inlimits.all(axis=1)

    elif method == 'river':
        rt = 'river'
        passed = (data_dict[rt] >= site_limits[0].loc[:, rt]) & (
        data_dict[rt] <= site_limits[0].loc[:, rt])  # todo check

    # todo sort out weighted options
    else:
        raise ValueError('should not get here')

    passed = pd.Series(index=nsmc_nums, data=passed)
    return passed
