# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/12/2017 9:56 AM
"""

from __future__ import division
from time import time
import numpy as np
import pandas as pd
import os
import multiprocessing
import logging
import psutil
import datetime
import socket
from traceback import format_exc
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.setup_reverse_modpath import setup_run_backward_modpath
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.extract_data import extract_back_data, save_forward_data
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.setup_forward_modpath import setup_run_forward_modpath
from waimak_extended_boundry.model_run_tools.model_setup.modpath_wrapper import get_cbc, get_cbc_mp
from waimak_extended_boundry.model_run_tools.metadata_managment.convergance_check import modpath_converged
import gc
import env


def define_source_from_forward(emulator_path, bd_type_path, indexes, return_packed_bits=False):
    """
    defines the source area array for a given integer array, values of the sourcea area array are counts and includes
    all particles which pass through cells flagged as an area of interest.
    this just extracts data from pre-run forward models
    :param emulator_path: path to the emulator (hdf) of forward modpath data
    :param bd_type: the boundary type assement from defineing particles, should be saved with the model run as a text
                    array.
    :param index: a dictionary of boolean arrays of areas of interest False delneates no interest, keys strings of any length
    :param return_packed_bits: bool if True return a boolean array as a packed bits array
    :return: dictionary of arrays of shape (rows, cols) with a particle count of the number of particles that
             originated in that cell that pass through the True index cell
    """
    # run some checks on inputs
    bd_type = np.loadtxt(bd_type_path).astype(int)
    assert isinstance(indexes, dict), 'indexes must be a dictionary'
    for key, idx in indexes.items():
        assert isinstance(idx, np.ndarray), 'index for {} must be a nd array'.format(key)
        assert idx.shape == (smt.layers, smt.rows, smt.cols), 'index for {} must be 3d'.format(key)
        assert idx.dtype == bool, 'index for {} must be some sort of boolean array'.format(key)

    # load emulator and initialize outdata
    print('loading emulator')
    t = time()
    emulator = pd.read_hdf(emulator_path)  # this keeps the structure of everything

    print('took {} s to load emulator'.format(time() - t))
    t = time()

    # get general area of interest
    print('calculating area of interest') # todo this could be improved by querrying the hdf store
    index = smt.get_empty_model_grid(True)
    for value in indexes.values():
        index += value

    idx = index > 0
    layers, rows, cols = np.meshgrid(range(smt.layers), range(smt.rows), range(smt.cols), indexing='ij')
    ids = ['{:02d}_{:03d}_{:03d}'.format(k, i, j) for k, i, j in zip(layers[idx], rows[idx], cols[idx])]
    temp = np.in1d(emulator.index.values, ids)
    emulator = emulator.loc[temp]
    gc.collect()
    print('took {} s to general identify area'.format(time() - t))

    # calculate source percentage
    all_outdata = {}
    for g, idx in indexes.items():
        ids = ['{:02d}_{:03d}_{:03d}'.format(k, i, j) for k, i, j in zip(layers[idx], rows[idx], cols[idx])]
        temp = np.in1d(emulator.index.values, ids)
        temp_emulator = emulator.loc[temp]
        temp_emulator = temp_emulator.reset_index()
        temp = temp_emulator.groupby('Particle_Group').aggregate({'fraction': np.sum})
        # I assume that I set particle group to the flattened cell number with only active cells hence the need
        # for the boundary type array, clever.

        # populate array
        outdata = smt.get_empty_model_grid(False)
        idx = bd_type.flatten() != -1
        outdata = outdata.flatten()
        temp_array = outdata[idx]
        temp_array[temp.index.values - 1] = temp.fraction.values
        outdata[idx] = temp_array
        outdata = outdata.reshape((smt.rows, smt.cols))
        if return_packed_bits:
            outdata = np.packbits(outdata > 0)
        all_outdata[g] = outdata

    return all_outdata


def define_source_from_backward(indexes, mp_ws, mp_name, cbc_file, root3_num_part=1, capt_weak_s=False, recalc=False,
                                return_packed_bits=False):
    """
    define the source area for an integer index, this creates and runs the models as well as extracts data
    :param indexes: a dictionary of boolean arrays, keys are strs of any length
    :param mp_ws: path for the modpath model
    :param mp_name: name of the modpath model
    :param cbc_file: cbc_file for the base modflow model of interest
    :param root3_num_part: the cubic root of the number of particles (placed evenly in the cell) e.g.
                           root3_num_part of 2 places 8 particles in each cell
    :param capt_weak_s: bool if True terminate particles at weak sources
    :param recalc: bool if True rerun the model even if it exists
    :param return_packed_bits: bool if True retun the data as packed boolean arrays
    :return:
    """
    assert isinstance(indexes, dict), 'indexes must be a dictionary'
    for key, idx in indexes.items():
        assert isinstance(idx, np.ndarray), 'index for {} must be a nd array'.format(key)
        assert idx.shape == (smt.layers, smt.rows, smt.cols), 'index for {} must be 3d'.format(key)
        assert idx.dtype == bool, 'index for {} must be some sort of integer array'.format(key)

    path_path = os.path.join(mp_ws, '{}.mppth'.format(mp_name))
    if not os.path.exists(path_path) or recalc:
        print('creating and running modpath model')
        # set up and run model
        setup_run_backward_modpath(mp_ws, mp_name, cbc_file, indexes,
                                   root3_num_part=root3_num_part, capt_weak_s=capt_weak_s)

    mapper_path = os.path.join(mp_ws, '{}_group_mapper.csv'.format(mp_name))
    outdata = extract_back_data(path_path, mapper_path, cbc_file.replace('.cbc', '.hds'),
                                return_packed_bits=return_packed_bits)
    return outdata


def _run_forward_em_one_mp(kwargs):
    """
    a wrapper to run the modpath as a multiprocessing
    (not that the multiprocessing was used as it was really IO and memory limited)
    used internally to run_forward_emulators
    :param kwargs:'model_id':
                  'mp_runs_dir':
                  'emulator_dir':
                  'modflow_dir':
                  'min_part':
                  'max_part':
                  'capt_weak_s':
                  'keep_org_files':
    :return:
    """
    model_id = kwargs['model_id']
    try:
        needed_keys = ['model_id', 'mp_runs_dir', 'emulator_dir', 'modflow_dir', 'min_part', 'max_part', 'capt_weak_s',
                       'keep_org_files']
        assert np.in1d(needed_keys, kwargs.keys()).all(), 'missing keys {}'.format(
            set(needed_keys) - set(kwargs.keys()))
        mp_ws = os.path.join(kwargs['mp_runs_dir'], model_id)
        outpath = os.path.join(kwargs['emulator_dir'], model_id + '.hdf')
        mp_name = '{}_forward'.format(model_id)
        cbc_path = get_cbc(model_id, kwargs['modflow_dir'])
        setup_run_forward_modpath(cbc_path, mp_ws, mp_name,
                                  min_part=kwargs['min_part'], max_part=kwargs['max_part'],
                                  capt_weak_s=kwargs['capt_weak_s'])
        path_path = os.path.join(mp_ws, '{}.mppth'.format(mp_name))
        save_forward_data(path_path, outpath)

        if not kwargs['keep_org_files']:
            # delete large files to save memory
            for end in ['.mpend', '.mppth', '.loc', '.mpbas']:  # all others are either needed or really tiny
                os.remove(os.path.join(mp_ws, '{}{}'.format(mp_name, end)))
        success = 'converge' if modpath_converged(
            os.path.join(mp_ws, '{}.mplst'.format(mp_name))) else 'did not converge'
    except Exception as val:
        mp_name = '{}_forward'.format(model_id)
        success = format_exc().replace('\n', '')

    return mp_name, success


def run_forward_emulators(model_ids, results_dir, modflow_dir, keep_org_files=True, min_part=1, max_part=None,
                          capt_weak_s=False, notes=''):
    """
    runs the forward emulators for the model ids
    :param model_ids: a list of model ids to run 'NsmcReal{nsmc_num:06d}'
    :param results_dir: the dir to put the modpath models and the results
    :param modflow_dir: the dir to put the necissary modflow models for modpath
    :param keep_org_files: bool if True keep all modpath files else delete the big ones
    :param min_part: the minimum number of particles for each cell
    :param max_part: the maximum number of particles for each cell
    :param capt_weak_s: bool if True terminate particles at weak sink
    :param notes: str any notes to be saves as a readme
    :return:
    """
    # set up a function to run all the emulators for the forward runs
    model_ids = np.atleast_1d(model_ids)
    emulator_dir = os.path.join(results_dir, 'forward_data')
    mp_runs_dir = os.path.join(results_dir, 'forward_runs')

    for path in [modflow_dir, results_dir, emulator_dir, mp_runs_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    with open(os.path.join(results_dir, 'README.txt'), 'w') as f:
        f.write(notes)

    input_kwargs = []
    for model_id in model_ids:
        temp = {'model_id': model_id,
                'mp_runs_dir': mp_runs_dir,
                'emulator_dir': emulator_dir,
                'modflow_dir': modflow_dir,
                'min_part': min_part,
                'max_part': max_part,
                'capt_weak_s': capt_weak_s,
                'keep_org_files': keep_org_files}

        input_kwargs.append(temp)
    t = time()
    # multiprocess the running of things
    outputs = []
    for i, kwarg in enumerate(input_kwargs):
        print('starting {} of {}. weak_sink?, {}'.format(i + 1, len(input_kwargs), capt_weak_s))
        outputs.append(_run_forward_em_one_mp(kwarg))
    now = datetime.datetime.now()
    with open(
            "{}/forward_run_log/{}_forward_modpath_{:02d}_{:02d}_{:02d}_{:02d}.txt".format(env.log_dir, now.year, now.month,
                                                                                           now.day, now.hour,
                                                                                           now.minute), 'w') as f:
        f.write(str(notes) + '\n')
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in outputs]
        f.writelines(wr)
    with open("{}/metadata.txt".format(results_dir), 'w') as f:
        f.write(str(notes) + '\n')
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in outputs]
        f.writelines(wr)
    print('{} runs completed in {} minutes'.format(len(model_ids), ((time() - t) / 60)))


def get_all_cbcs(model_ids, modflow_dir, sleep_time=1, recalc=False):
    """
    a quick multiprocessing wrapper to run all the NSMC realisations I need
    note modifying the modpath model is not implmented, this only runs on the base model
    :param model_ids: list of model ids to pass to the get cbc
    :param modflow_dir: overall directory to save everything
    :param sleep_time: the sleep time between printing number of runs left (mostly for debugging purposes)
    :return:
    """
    input_kwargs = []
    for model_id in model_ids:
        input_kwargs.append({'model_id': model_id,
                             'base_dir': modflow_dir,
                             'recalc': recalc})

    multiprocessing.log_to_stderr(logging.DEBUG)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map_async(get_cbc_mp, input_kwargs)
    pool_outputs = results.get()
    pool.close()  # no more tasks
    pool.join()
    with open(os.path.join(modflow_dir, 'metadata.txt'), 'w') as f:
        wr = ['{}: {}\n'.format(e[0], e[1]) for e in pool_outputs]
        f.writelines(wr)


def start_process():
    print('Starting', multiprocessing.current_process().name, 'pid: {}'.format(os.getpid()))
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def get_modeflow_dir_for_source(version=1):
    """
    somethign to keep the modflow dirs for the modpath sernarios straight
    :param version: which version (so far only 1 version)
    :return:
    """
    if version == 1:
        raise NotImplementedError ('modflow dir for source is not implemented')
    else:
        raise NotImplementedError('version {} has not been set up'.format(version))
    return path


def get_base_results_dir(mode, comp):
    """
    keeps track fo the base dirs for the forward and backward particles
    :param mode:'forward' or 'backward'
    :param comp: the socket host name
    :return:
    """

    raise NotImplementedError('this is a hold out from previous modelling, it is not currently implmented')

    if mode == 'forward' and comp == 'GWATER02':
        out = r"D:\mh_waimak_models_from_c\modpath_forward_base"
    elif mode == 'backward' and comp == 'GWATER02':
        out = r"D:\mh_waimak_models_from_c\modpath_reverse_base"
    elif mode == 'cust' and comp == 'GWATER02':
        out = r'D:\mh_waimak_models_from_c\cust_data'
    elif mode == 'forward' and comp == 'RDSProd03':
        out = r"D:\mh_waimak_models\modpath_forward_base"
    elif mode == 'backward' and comp == 'RDSProd03':
        out = r"D:\mh_waimak_models\modpath_reverse_base"
    elif mode == 'cust' and comp == 'RDSProd03':
        out = r'D:\mh_waimak_models\cust_data'
    else:
        raise ValueError('unexpected (mode, comp): ({},{})'.format(mode, comp))
    return out


def get_forward_emulator_paths(model_ids, weak_sink=False):
    """
    get dictionary of forward emulator paths, this was a helper to keep track of particle tracking runs, these data
    were not stored and this is a bit historical.  it could be re-purposed for future work though
    :param model_ids: model_ids
    :param weak_sink: Bool if True capture at weak sink
    :return: {model_id: (emulator_path, bnd_type_path)}
    """
    model_ids = np.atleast_1d(model_ids)
    forward_base = get_base_results_dir(mode='forward', comp=socket.gethostname())
    if weak_sink:
        sink = 'weak_sinks'
    else:
        sink = 'strong_sinks'
    em_paths = [os.path.join(forward_base, sink, 'forward_data', '{}.hdf'.format(mid)) for mid in model_ids]
    bnd_types = [os.path.join(forward_base, sink, 'forward_runs', mid,
                              '{}_forward_bnd_type.txt'.format(mid)) for mid in model_ids]
    out = {mid: (em, bd) for mid, em, bd in zip(model_ids, em_paths, bnd_types)}
    return out


if __name__ == '__main__':
    # this looks good
    # check layer 7 in the area of no data
    import matplotlib.pyplot as plt

    test_type = 3
    if test_type == 1:
        temp_index = smt.shape_file_to_model_array(r"C:\Users\MattH\Downloads\test_area.shp", 'Id', True)
        index = smt.get_empty_model_grid(True).astype(bool)
        index[1] = np.isfinite(temp_index)
        index = smt.get_empty_model_grid(True).astype(bool)
        index = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\rough_chch.shp".format(smt.sdp), 'Id', True)[
            np.newaxis].repeat(11, axis=0)
        index = np.isfinite(index).astype(int)
        bd_type = (
        r"C:\mh_waimak_models\modpath_forward_base\strong_sinks\forward_runs\NsmcReal-00001.hdf\NsmcReal-00001_forward_bnd_type.txt")
        outdata = define_source_from_forward(
            r"C:\mh_waimak_models\modpath_forward_base\strong_sinks\forward_data\NsmcReal-00001", bd_type,
            index.astype(int))
        outdata[1][outdata[1] == 0] = np.nan
        smt.plt_matrix(np.log10(outdata[1]), base_map=True)
        plt.show()

    elif test_type == 2:
        import pickle

        index = smt.get_empty_model_grid(True).astype(bool)
        index = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\rough_chch.shp".format(smt.sdp), 'Id', True)[
            np.newaxis]
        index2 = np.isfinite(index).repeat(11, axis=0)
        index2[1:, :, :] = False
        index1 = np.full((smt.layers, smt.rows, smt.cols), False)
        index1[6] = index
        indexes = {'layer0_chch': index2, 'layer7_chch': index1}
        outdata = define_source_from_backward(indexes,
                                              r"C:\Users\MattH\Downloads\test_back",
                                              'test_back',
                                              get_cbc('NsmcBase', get_modeflow_dir_for_source()),
                                              recalc=False)
        pickle.dump(outdata, open(r"C:\Users\MattH\Downloads\testback_zones.p", 'w'))
        print('done')

    if test_type==3:
        nsmc_nums = [2971, 3116, 3310, 3378, 3456, 3513, 3762, 3910]
        model_ids = ['NsmcReal{:06d}'.format(e) for e in nsmc_nums]
        get_all_cbcs(model_ids, get_modeflow_dir_for_source(), recalc=True)
