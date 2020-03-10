"""
Author: matth
Date Created: 31/05/2017 9:19 AM
"""

from __future__ import division
from users.MH.Waimak_modeling.supporting_data_path import sdp
import flopy
import os
from stream_depletion_model_setup import setup_and_run_stream_dep
import transient_inputs as tinputs
import numpy as np


def get_starting_heads_for_sd(model_version='m_strong_vert', recalc=False):
    """
    return starting heads
    :return:
    """
    heads_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-starting_hds_run/{v}-starting_hds_run.hds".format(
        p=sdp, v=model_version)
    if os.path.exists(heads_path) and not recalc:
        hds = flopy.utils.HeadFile(heads_path)
        out_heads = hds.get_data(hds.get_kstpkper()[-1])  # use the last timestep
        return out_heads

    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    spv = {'nper': 4,  #
           'perlen': 30,
           'nstp': 2,
           'steady': [True, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: None,
                       1: 7,
                       2: 8,
                       3: 9,
                       }

    base_path = "{}/base_model_runs/base_stream_depletion_mf".format(sdp)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    kwargs = {'name': 'starting_hds_run',
              'base_dir': base_path,
              'stress_vals': spv,
              'stress_to_month': stress_to_month,
              'wells_to_turn_on': {0: []},
              'ss': ss,
              'sy': sy,
              'silent': False,
              'solver': 'pcg',
              'model_version': model_version,
              'model': 'mfnwt'}  # solver will be changed internally to nwt

    setup_and_run_stream_dep(**kwargs)
    hds = flopy.utils.HeadFile(heads_path)
    out_heads = hds.get_data(hds.get_kstpkper()[-1])  # use the last timestep
    return out_heads


def get_fully_nat_str_dep150_base_path(model_version='m_strong_vert', recalc=False, add_h20_to_cust=False):
    """
    function to return the model path for the fully naturalized fully transient run and ensure that the model has been
    run this run has a ss and then full year run from july to end of june  ## for SD 150
    :param recalc:
    :return: path to the model (without and extension)
    """

    if add_h20_to_cust:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd150_1cu_cust/{v}-fully_naturalized_sd150_1cu_cust".format(
            p=sdp, v=model_version)
    else:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd150_run/{v}-fully_naturalized_sd150_run".format(
            p=sdp, v=model_version)
    if os.path.exists('{}.hds'.format(model_path)) and not recalc:
        return model_path

    # set up fully transient run
    base_path = "{}/base_model_runs/base_stream_depletion_mf".format(sdp)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    start_hds = get_starting_heads_for_sd()
    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    spv = {'nper': 5,
           'perlen': 30,
           'nstp': 2,
           'steady': [False, False, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: 10,
                       1: 11,
                       2: 12,
                       3: 1,
                       4: 2}
    if add_h20_to_cust:
        name = 'fully_naturalized_sd150_1cu_cust'
        add_h20 = 86400
    else:
        name = 'fully_naturalized_sd150_run'
        add_h20 = 0
    kwargs = {'name': name,
              'base_dir': base_path,
              'stress_vals': spv,
              'stress_to_month': stress_to_month,
              'wells_to_turn_on': {0: []},
              'ss': ss,
              'sy': sy,
              'silent': False,
              'solver': 'pcg',
              'model': 'mfnwt',  # solver will be changed internally to nwt
              'start_heads': start_hds,
              'flow_to_add_to_cust': add_h20,
              'model_version': model_version}

    setup_and_run_stream_dep(**kwargs)

    return model_path


def get_fully_nat_str_dep7_base_path(model_version='m_strong_vert', recalc=False, add_h20_to_cust=False):
    """
    function to return the model path for the fully naturalized fully transient run and ensure that the model has been
    run this run has a ss and then full year run from july to end of june  ## for SD 150
    :param recalc:
    :return: path to the model (without and extension)
    """
    if add_h20_to_cust:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd7_1cu_cust/{v}-fully_naturalized_sd7_1cu_cust".format(
            p=sdp, v=model_version)
    else:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd7_run/{v}-fully_naturalized_sd7_run".format(
            p=sdp, v=model_version)
    if os.path.exists('{}.hds'.format(model_path)) and not recalc:
        return model_path

    # set up fully transient run
    base_path = "{}/base_model_runs/base_stream_depletion_mf".format(sdp)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    spv = {'nper': 8,
           'perlen': 1,
           'nstp': 1,
           'steady': [True, False, False, False, False, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: None,
                       1: None,
                       2: None,
                       3: None,
                       4: None,
                       5: None,
                       6: None,
                       7: None}

    if add_h20_to_cust:
        name = 'fully_naturalized_sd7_1cu_cust'
        add_h20 = 86400
    else:
        name = 'fully_naturalized_sd7_run'
        add_h20 = 0
    kwargs = {'name': name,
              'base_dir': base_path,
              'stress_vals': spv,
              'stress_to_month': stress_to_month,
              'wells_to_turn_on': {0: []},
              'ss': ss,
              'sy': sy,
              'silent': False,
              'solver': 'pcg',
              'flow_to_add_to_cust': add_h20,
              'model': 'mfnwt',  # solver will be changed internally to nwt
              'model_version': model_version}

    setup_and_run_stream_dep(**kwargs)

    return model_path


def get_fully_nat_str_dep30_base_path(model_version='m_strong_vert', recalc=False, add_h20_to_cust=False):
    """
    function to return the model path for the fully naturalized fully transient run and ensure that the model has been
    run this run has a ss and then full year run from july to end of june  ## for SD 150
    :param recalc:
    :return: path to the model (without and extension)
    """
    if add_h20_to_cust:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd30_1cu_cust/{v}-fully_naturalized_sd30_1cu_cust".format(
            p=sdp, v=model_version)
    else:
        model_path = "{p}/base_model_runs/base_stream_depletion_mf/{v}-fully_naturalized_sd30_run/{v}-fully_naturalized_sd30_run".format(
            p=sdp, v=model_version)
    if os.path.exists('{}.hds'.format(model_path)) and not recalc:
        return model_path

    # set up fully transient run
    base_path = "{}/base_model_runs/base_stream_depletion_mf".format(sdp)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    log_dir = '{}/logging'.format(base_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    ss = np.percentile(tinputs.get_ss_dist(), 50)
    sy = np.percentile(tinputs.get_sy_dist(), 50)

    spv = {'nper': 11,
           'perlen': 3,
           'nstp': 1,
           'steady': [True, False, False, False, False, False, False, False, False, False, False],
           'tsmult': 1.1}

    stress_to_month = {0: None,
                       1: None,
                       2: None,
                       3: None,
                       4: None,
                       5: None,
                       6: None,
                       7: None,
                       8: None,
                       9: None,
                       10: None}

    if add_h20_to_cust:
        name = 'fully_naturalized_sd30_1cu_cust'
        add_h20 = 86400
    else:
        name = 'fully_naturalized_sd30_run'
        add_h20 = 0
    kwargs = {'name': name,
              'base_dir': base_path,
              'stress_vals': spv,
              'stress_to_month': stress_to_month,
              'wells_to_turn_on': {0: []},
              'ss': ss,
              'sy': sy,
              'silent': False,
              'solver': 'pcg',
              'flow_to_add_to_cust': add_h20,
              'model': 'mfnwt',  # solver will be changed internally to nwt
              'model_version': model_version}

    setup_and_run_stream_dep(**kwargs)

    return model_path


if __name__ == '__main__':
    # get_starting_heads_for_sd(recalc=True)
    get_fully_nat_str_dep150_base_path(recalc=True, add_h20_to_cust=True)
    get_fully_nat_str_dep7_base_path(recalc=True, add_h20_to_cust=True)
    get_fully_nat_str_dep30_base_path(recalc=True, add_h20_to_cust=True)
