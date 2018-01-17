# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/01/2018 2:21 PM
"""

from __future__ import division
from core import env
import numpy as np
import netCDF4 as nc
import os
import socket
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.source_delmination.source_raster import \
    define_source_from_backward, define_source_from_forward, get_modeflow_dir_for_source, get_base_results_dir, \
    get_cbc, get_forward_emulator_paths


# set up all single zone delination maps (e.g. streams/ water supply wells)
# todo change indexes to be a dictionary of indexes to avoid overlaps?
# todo how to handle the cust/+- ashley problem
def create_single_zone_indexs():
    """
    create an index for all zones
    :return: {id: index}
    """
    raise NotImplementedError


def create_zones(model_ids, outpath, root_num_part, recalc=False, recalc_backward_tracking=False):
    """
    set up something to make a dictionary of mappers and return out put, make it save and/or load it so that I don't
    need to keep re-calculating results as netcdf (with geocoding) do both strong/weak forward/backward
    :return:
    """
    if recalc_backward_tracking:
        recalc = True
    model_ids = np.atleast_1d(model_ids)

    if os.path.exists(outpath) and not recalc:
        outdata = nc.Dataset(outpath)
        return outdata

    indexes = create_single_zone_indexs()
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = get_base_results_dir('backward', socket.gethostname())

    # backward weak
    back_weaks = []
    for mid in model_ids:
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'weak', mid),
                                           mp_name='{}_weak'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=True,
                                           recalc=recalc_backward_tracking)
        back_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_weak_back = _amalg_backward(back_weaks, root_num_part, indexes)

    # backward strong
    back_strongs = []
    for mid in model_ids:
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'strong', mid),
                                           mp_name='{}_strong'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=False,
                                           recalc=recalc_backward_tracking)
        back_strongs.append(temp)
    # amalgamate data from realsiations mean and sd?

    # forward weak
    forward_weaks = []
    forward_weaks_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, True)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type=path[1], indexes=indexes)
        temp2 = np.loadtxt(path[1].replace('_forward_bnd_type.txt','_num_parts.txt')) # load number of particles
        forward_weaks_num_parts.append(temp2)
        forward_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?

    # forward strong
    forward_strongs = []
    forward_strongs_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type=path[1], indexes=indexes)
        forward_strongs.append(temp)
        temp2 = np.loadtxt(path[1].replace('_forward_bnd_type.txt','_num_parts.txt')) # load number of particles
        forward_strongs_num_parts.append(temp2)
    # amalgamate data from realsiations

    # save the data

    # re-open the data and return
    raise NotImplementedError


def plot_sources():
    """
    take the output of create zones and plot it up
    :return:
    """

    raise NotImplementedError


def _save_source_nc(forward_weak, forward_strong, backward_strong, backward_weak, model_ids):
    """
    save the thing as a netcdf file
    dimensions for data(direction, source behavior, row, col) variables of location ids
    :return:
    """
    # todo add the number of particles for forward model? and the root num particles for the backward
    raise NotImplementedError


def _amalg_forward(): #todo

    raise NotImplementedError

def _amalg_backward(data, root_num_part, indexes):
    """

    :param data: the list of dictionaries for the data
    :param root_num_part: the root number of particles generated
    :param indexes: the indexes going forward
    :return: {id: map} id should match amalg forward
    """
    out = {}
    for name in indexes.keys():
        temp = [e[name] for e in data] #todo normalise by number of particles and figure out nan behavior, which I don't think should be a problem
        out['{}_u'.format(name)] = np.mean(temp, axis=0)
        out['{}_sd'.format(name)] = np.std(temp, axis=0)
        out['{}_sum'.format(name)] = np.sum(temp, axis=0)
        out['{}_min'.format(name)] = np.min(temp, axis=0)
        out['{}_max'.format(name)] = np.min(temp, axis=0)
        out['{}_max'.format(name)] = np.min(temp, axis=0) #todo add percentiles anything else? start here



    raise NotImplementedError
