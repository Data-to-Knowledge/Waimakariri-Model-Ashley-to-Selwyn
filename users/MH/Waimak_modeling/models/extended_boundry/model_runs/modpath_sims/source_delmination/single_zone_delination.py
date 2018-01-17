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
import sys
import datetime
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.source_delmination.source_raster import \
    define_source_from_backward, define_source_from_forward, get_modeflow_dir_for_source, get_base_results_dir, \
    get_cbc, get_forward_emulator_paths


# set up all single zone delination maps (e.g. streams/ water supply wells)
# todo how to handle the cust/+- ashley problem? manually? that would be a pain in the ass I;m not sure what to do here...
# I could run one for all of the cust and add the particles... would mess up my numbering normalisation...

def create_single_zone_indexs(): #todo
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
    amalg_strong_back = _amalg_backward(back_strongs, root_num_part, indexes)

    # forward weak
    forward_weaks = []
    forward_weaks_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, True)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type=path[1], indexes=indexes)
        temp2 = np.loadtxt(path[1].replace('_forward_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_weaks_num_parts.append(temp2)
        forward_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_weak_forward = _amalg_forward(forward_weaks, forward_weaks_num_parts, indexes)

    # forward strong
    forward_strongs = []
    forward_strongs_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type=path[1], indexes=indexes)
        forward_strongs.append(temp)
        temp2 = np.loadtxt(path[1].replace('_forward_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_strongs_num_parts.append(temp2)
    # amalgamate data from realsiations
    amalg_strong_forward = _amalg_forward(forward_strongs, forward_strongs_num_parts, indexes)

    # save the data
    outdata = _save_source_nc(outpath=outpath,
                              forward_weak=amalg_weak_forward,
                              forward_strong=amalg_strong_forward,
                              backward_strong=amalg_strong_back,
                              backward_weak=amalg_weak_back,
                              model_ids=model_ids,
                              root_num_part=root_num_part,
                              forward_part_nums=forward_strongs_num_parts)
    return outdata


def plot_sources():  # todo play in arc first
    """
    take the output of create zones and plot it up
    :return:
    """

    raise NotImplementedError


def _save_source_nc(outpath, forward_weak, forward_strong, backward_strong, backward_weak, model_ids, root_num_part,
                    forward_part_nums):
    """
    save the thing as a netcdf file
    dimensions for data(direction, source behavior, amalg_type, row, col) variables of location ids
    :return:
    """
    # some assertions
    assert all([isinstance(e, dict) for e in
                [forward_strong, forward_weak, backward_strong, backward_weak]]), 'all data must be dict'
    temp = set(forward_weak.keys()) == set(forward_strong.keys()) == set(backward_weak.keys()) == set(
        backward_strong.keys())
    assert temp, 'all inputs must have the same keys'

    amalg_types = list(set([e.split('~')[-1] for e in forward_weak.keys()]))
    sites = list(set([e.split('~')[0] for e in forward_weak.keys()]))

    outfile = nc.Dataset(outpath, 'w')
    x, y = smt.get_model_x_y(False)
    # create dimensions
    outfile.createDimension('latitude', len(y))
    outfile.createDimension('longitude', len(x))
    outfile.createDimension('sim_dir', 2)
    outfile.createDimension('capt_behav', 2)
    outfile.createDimension('amalg_type', len(amalg_types))

    # create variables

    proj = outfile.createVariable('crs', 'i1')  # this works really well...
    proj.setncatts({'grid_mapping_name': "transverse_mercator",
                    'scale_factor_at_central_meridian': 0.9996,
                    'longitude_of_central_meridian': 173.0,
                    'latitude_of_projection_origin': 0.0,
                    'false_easting': 1600000,
                    'false_northing': 10000000,
                    })

    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'latitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_y_coordinate'})
    lat[:] = y

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lon.setncatts({'units': 'NZTM',
                   'long_name': 'longitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_x_coordinate'})
    lon[:] = x

    sim_dir = outfile.createVariable('sim_dir', str, ('sim_dir',))
    sim_dir.setncatts({'long_name': 'simulation_direction'})
    sim_dir[:] = ['forward', 'backward']

    capt_bev = outfile.createVariable('capt_behav', str, ('capt_behav',))
    capt_bev.setncatts({'long_name': 'particle_capture_behavior'})
    capt_bev[:] = ['strong', 'weak']

    amalg_type = outfile.createVariable('amalg_type', str, ('amalg_type',))
    amalg_type.setncatts({'long_name': 'realisation_amalgamation_types'})
    amalg_type[:] = amalg_types

    # add the number of particles for forward model (same amalgimation types)
    temp_data = {}
    _add_data_variations(temp_data, forward_part_nums, 'forward_part_nums')
    temp_var = outfile.createVariable('forward_part_nums', int,
                                      ('amalg_type', 'latitude', 'longitude'),
                                      fill_value=-1)
    temp_var.setncatts({'units': 'none',
                        'long_name': 'number of particles released for forward run',
                        'missing_value': -1})

    temp_data = np.full((len(amalg_types), smt.rows, smt.cols), np.nan)
    # some for loop trickery
    for a, at in enumerate(amalg_types):
        temp_var[a, :, :] = forward_strong['forward_part_nums~{}'.format(at)]
    temp_var[:] = temp_data

    # location add the data
    for site in sites:
        temp_var = outfile.createVariable(site, float,
                                          ('sim_dir', 'capt_behav', 'amalg_type', 'latitude', 'longitude'),
                                          fill_value=np.nan)
        temp_var.setncatts({'units': 'fraction',
                            'long_name': site,
                            'missing_value': np.nan,
                            'comments': 'fraction of particles released in cell (forward) or '
                                        'fraction of particles released for site (backward)'})

        temp_data = np.full((2, 2, len(amalg_types), smt.rows, smt.cols), np.nan)
        # some for loop trickery
        for a, at in enumerate(amalg_types):
            temp_var[0, 0, a, :, :] = forward_strong['{}~{}'.format(site, at)]
            temp_var[0, 1, a, :, :] = forward_weak['{}~{}'.format(site, at)]
            temp_var[1, 0, a, :, :] = backward_strong['{}~{}'.format(site, at)]
            temp_var[1, 1, a, :, :] = backward_weak['{}~{}'.format(site, at)]
        temp_var[:] = temp_data

    outfile.description = ('source zones for single sources')
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'script: {}'.format(sys.argv[0])
    outfile.backward_num_parts = root_num_part ** 3
    outfile.model_ids = model_ids  # todo check that it works to assign a list to a description
    outfile.close()

    outfile = nc.Dataset(outpath)  # reload to avoid write behavior
    return outfile


def _amalg_forward(data, num_parts, indexes):
    out = {}
    for name in indexes.keys():
        temp = [e[name] / np for e, np in
                zip(data, num_parts)]  # todo figure out nan behavior, which I don't think should be a problem
        _add_data_variations(out, temp, name)

    return out


def _amalg_backward(data, root_num_part, indexes):
    """

    :param data: the list of dictionaries for the data
    :param root_num_part: the root number of particles generated
    :param indexes: the indexes going forward
    :return: {id: map} id should match amalg forward
    """
    out = {}
    for name in indexes.keys():
        temp = [e[name] / root_num_part ** 3 for e in
                data]  # todo figure out nan behavior, which I don't think should be a problem
        _add_data_variations(out, temp, name)

    return out


def _add_data_variations(out, temp, name):
    """
    quick wrapper to ensure the same behaviour between forward and backward amalgamations
    :param out:
    :param temp:
    :param name:
    :return:
    """
    out['{}~u'.format(name)] = np.mean(temp, axis=0)
    out['{}~sd'.format(name)] = np.std(temp, axis=0)
    out['{}~sum'.format(name)] = np.sum(temp, axis=0)
    out['{}~min'.format(name)] = np.min(temp, axis=0)
    out['{}~max'.format(name)] = np.min(temp, axis=0)
    out['{}~5th'.format(name)] = np.percentile(temp, 5, axis=0)
    out['{}~25th'.format(name)] = np.percentile(temp, 25, axis=0)
    out['{}~50th'.format(name)] = np.percentile(temp, 50, axis=0)
    out['{}~75th'.format(name)] = np.percentile(temp, 75, axis=0)
    out['{}~95th'.format(name)] = np.percentile(temp, 95, axis=0)
    temp2 = [e > 0 for e in temp]
    out['{}~all'.format(name)] = np.all(temp2, axis=0)
    out['{}~any'.format(name)] = np.any(temp2, axis=0)
