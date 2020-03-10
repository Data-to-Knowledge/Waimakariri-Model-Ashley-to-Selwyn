"""
Author: matth
Date Created: 1/06/2017 11:58 AM
"""

from __future__ import division
import flopy
import pandas as pd
from base_str_dep_runs import get_fully_nat_str_dep150_base_path, get_fully_nat_str_dep7_base_path, \
    get_fully_nat_str_dep30_base_path
import users.MH.Waimak_modeling.model_tools as mt
import socket
import glob
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, results_dir, base_mod_dir2
import itertools
import numpy as np
from core.ecan_io import rd_sql, sql_db
import multiprocessing
import logging


def calc_stream_dep(model_path, mode='sd150', add_h20_cust=False):
    """
    calculate a dataframe for the stream depletion
    :param model_path:
    :param mode: either sd150 or sd7,  this defines which base mod path to use.
    :return:
    """

    if mode == 'sd150':
        baseline_path = get_fully_nat_str_dep150_base_path(add_h20_to_cust=add_h20_cust)
    elif mode == 'sd7':
        baseline_path = get_fully_nat_str_dep7_base_path(add_h20_to_cust=add_h20_cust)
    elif mode == 'sd30':
        baseline_path = get_fully_nat_str_dep30_base_path(add_h20_to_cust=add_h20_cust)
    else:
        raise ValueError('unexpected stream depletion mode {}'.format(mode))

    str_points = mt.get_str_sample_points()
    idx = ['sd' in e for e in str_points]
    str_points = list(np.array(str_points)[np.array(idx)])

    drn_points = mt.get_drn_samp_pts()
    idx = ['swaz' in e for e in drn_points]
    drn_points = list(np.array(drn_points)[np.array(idx)])
    drn_points.remove('ashley_swaz')
    drn_points.remove('num7drain_swaz')
    drn_points.remove('str_eyre_swaz')

    if mode == 'sd150':
        kstpkpers_base = list(itertools.product(range(2), range(0, 5)))
        kstpkpers_sd = list(itertools.product(range(2), range(0, 5)))
        integrater = 15
    elif mode == 'sd7':
        kstpkpers_base = list(itertools.product(range(1), range(1, 8)))
        kstpkpers_sd = list(itertools.product(range(1), range(1, 8)))
        integrater = 1
    elif mode == 'sd30':
        kstpkpers_base = list(itertools.product(range(1), range(1, 11)))
        kstpkpers_sd = list(itertools.product(range(1), range(1, 11)))
        integrater = 3
    else:
        raise ValueError('unexpected stream depletion mode {}'.format(mode))

    str_base = mt.streamflow_for_kskps(baseline_path,
                                       kstpkpers=kstpkpers_base,
                                       drn_points=drn_points, str_points=str_points)
    str_sd = mt.streamflow_for_kskps(model_path, kstpkpers=kstpkpers_sd,
                                     drn_points=drn_points, str_points=str_points)

    outdata = (str_base.sum() - str_sd.sum()) * integrater

    # cumulative well abstraction budget
    budget = flopy.utils.MfListBudget('{}.list'.format(model_path))
    temp = pd.DataFrame(budget.get_cumulative('WELLS_OUT'))
    abs_vol = temp[(temp['stress_period'] == temp['stress_period'].max()) &
                   (temp['time_step'] == temp['time_step'].max())]['WELLS_OUT'].iloc[0]
    # convert from l3 to m3

    if abs_vol == 0:
        outdata *= np.nan
    else:
        outdata *= 1 / abs_vol
    return outdata


def calc_str_dep_all_wells(base_path, mode='sd150', add_h20_cust=False):
    all_paths = glob.glob('{}/*/*.nam'.format(base_path))
    all_paths = [e.strip('.nam') for e in all_paths]
    wells = ['{}/{}'.format(e.split('_')[-2], e.split('_')[-1]) for e in all_paths]
    outdata = {}
    for well, path in zip(wells, all_paths):
        outdata[well] = calc_stream_dep(path, mode=mode, add_h20_cust=add_h20_cust)
    outdata = pd.DataFrame(outdata).transpose()

    return outdata


def join_sd_results(version='sd150',add_water_to_cust=False):
    """
    a script to join up all the stream depletion results.  it must be first run on the remote machines to get their
    results to the server and then on my machine to compile the results
    :param version: either 'sd150', 'sd7', or 'all'
    :return:
    """
    _comp = socket.gethostname()

    if add_water_to_cust:
        datasd7_path = '{}/stream_depletion/stream_depletion_7/stream_d7_data_1cu_cust_'.format(results_dir)
        datasd150_path = '{}/stream_depletion/stream_depletion_150/stream_d150_data_1cu_cust_'.format(results_dir)
        datasd30_path = '{}/stream_depletion/stream_depletion_30/stream_d30_data_1cu_cust_'.format(results_dir)
        root_path = "{}/well_by_well_str_dep".format(base_mod_dir2)
        extra_path = '_1cu_added_to_cust'
    else:
        datasd7_path = '{}/stream_depletion/stream_depletion_7/stream_d7_data_'.format(results_dir)
        datasd150_path = '{}/stream_depletion/stream_depletion_150/stream_d150_data_'.format(results_dir)
        datasd30_path = '{}/stream_depletion/stream_depletion_30/stream_d30_data_'.format(results_dir)
        root_path = '{}/well_by_well_str_dep'.format(base_mod_dir)
        extra_path = ''

    if _comp == 'DHI-Runs02' or _comp == 'GWATER02':
        if version == 'sd150' or version == 'all':
            temp_path = '{}_sd150{}'.format(root_path,extra_path)
            datasd150 = calc_str_dep_all_wells(temp_path, mode='sd150', add_h20_cust=add_water_to_cust)
            datasd150.to_csv('{}{}.csv'.format(datasd150_path, _comp))
        if version == 'sd7' or version == 'all':
            temp_path = '{}_sd7{}'.format(root_path,extra_path)
            datasd7 = calc_str_dep_all_wells(temp_path, mode='sd7', add_h20_cust=add_water_to_cust)
            datasd7.to_csv('{}{}.csv'.format(datasd7_path, _comp))
        if version == 'sd30' or version == 'all':
            temp_path = '{}_sd30{}'.format(root_path,extra_path)
            datasd30 = calc_str_dep_all_wells(temp_path, mode='sd30', add_h20_cust=add_water_to_cust)
            datasd30.to_csv('{}{}.csv'.format(datasd30_path, _comp))


    elif _comp == 'HP1639':
        # combine the two datasets and add other info (depth, CRC number, x, y )
        org_data = mt.get_original_well_data()
        well_details = rd_sql(**sql_db.wells_db.well_details)
        well_details = well_details.set_index('WELL_NO')
        well_details = well_details.loc[:, ['NZTMX', 'NZTMY']]
        well_details = well_details.rename(columns={'NZTMX': 'wells_x', 'NZTMY': 'wells_y'})

        if version == 'sd150' or version == 'all':
            flux_data = pd.DataFrame(mt.get_model_well_full_consented()['flux'])
            flux_data *= 12 / 5
            data_gwruns = pd.read_csv('{}GWATER02.csv'.format(datasd150_path))
            data_dhiruns = pd.read_csv('{}DHI-Runs02.csv'.format(datasd150_path))
            out_data = pd.concat((data_dhiruns, data_gwruns))
            out_data = out_data.set_index(out_data.keys()[0])

            # add original data
            out_data = pd.merge(out_data, org_data, how='left', left_index=True, right_index=True)
            out_data = out_data.drop(['flux', 'layer', 'row'], axis=1)
            out_data = pd.merge(out_data, well_details, how='left', left_index=True, right_index=True)
            out_data['mis_location'] = ((out_data['x'] - out_data['wells_x']) ** 2 + (
            out_data['y'] - out_data['wells_y']) ** 2) ** 0.5
            out_data = pd.merge(out_data, flux_data, how='left', left_index=True, right_index=True)

            out_data.to_csv('{}joined.csv'.format(datasd150_path))

        if version == 'sd7' or version == 'all':
            flux_data = pd.DataFrame(mt.get_model_well_max_rate()['flux'])
            data_gwruns = pd.read_csv('{}GWATER02.csv'.format(datasd7_path))
            data_dhiruns = pd.read_csv('{}DHI-Runs02.csv'.format(datasd7_path))
            out_data = pd.concat((data_dhiruns, data_gwruns))
            out_data = out_data.set_index(out_data.keys()[0])

            # add original data
            out_data = pd.merge(out_data, org_data, how='left', left_index=True, right_index=True)
            out_data = out_data.drop(['flux', 'layer', 'row'], axis=1)
            out_data = pd.merge(out_data, well_details, how='left', left_index=True, right_index=True)
            out_data['mis_location'] = ((out_data['x'] - out_data['wells_x']) ** 2 +
                                        (out_data['y'] - out_data['wells_y']) ** 2) ** 0.5
            out_data = pd.merge(out_data, flux_data, how='left', left_index=True, right_index=True)
            out_data.to_csv('{}joined.csv'.format(datasd7_path))

        if version == 'sd30' or version == 'all':
            flux_data = pd.DataFrame(mt.get_model_well_max_rate()['flux'])
            data_gwruns = pd.read_csv('{}GWATER02.csv'.format(datasd30_path))
            data_dhiruns = pd.read_csv('{}DHI-Runs02.csv'.format(datasd30_path))
            out_data = pd.concat((data_dhiruns, data_gwruns))
            out_data = out_data.set_index(out_data.keys()[0])

            # add original data
            out_data = pd.merge(out_data, org_data, how='left', left_index=True, right_index=True)
            out_data = out_data.drop(['flux', 'layer', 'row'], axis=1)
            out_data = pd.merge(out_data, well_details, how='left', left_index=True, right_index=True)
            out_data['mis_location'] = ((out_data['x'] - out_data['wells_x']) ** 2 +
                                        (out_data['y'] - out_data['wells_y']) ** 2) ** 0.5
            out_data = pd.merge(out_data, flux_data, how='left', left_index=True, right_index=True)
            out_data.to_csv('{}joined.csv'.format(datasd30_path))

    else:
        raise ValueError('unexpected computer {}'.format(_comp))


if __name__ == '__main__':
    run_type = 2
    if run_type == 2:
        join_sd_results('all',True)
