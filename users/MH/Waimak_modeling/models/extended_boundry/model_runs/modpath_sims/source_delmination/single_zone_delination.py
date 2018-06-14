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
import flopy_mh as flopy
import pandas as pd
from copy import deepcopy
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.source_delmination.source_raster import \
    define_source_from_backward, define_source_from_forward, get_modeflow_dir_for_source, get_base_results_dir, \
    get_cbc, get_forward_emulator_paths
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_from_streams import \
    get_samp_points_df, _get_sw_samp_pts_dict, _get_flux_flow_arrays
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set


def create_single_zone_indexs():
    """
    create an indexs for all zones for the stream and WDC wells
    :return: {id: index}
    """
    indexes = {}
    # streams
    str_data = get_samp_points_df()
    str_data_use = str_data.loc[str_data.m_type == 'source']
    str_dict = _get_sw_samp_pts_dict()
    for idx in str_data_use.index:
        sfr_array, drn_array = _get_flux_flow_arrays(idx, str_dict, str_data)
        temp_out = smt.get_empty_model_grid(True).astype(bool)
        if sfr_array is not None:
            temp_out[0] = sfr_array.astype(bool) | temp_out[0]
        if drn_array is not None:
            temp_out[0] = drn_array.astype(bool) | temp_out[0]
        indexes[idx] = temp_out

    # WDC water supply wells (groups of individual wells)
    wdc_wells = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\WDC_wells.csv")
    all_wells = get_all_well_row_col()

    for z in set(wdc_wells.Zone):
        temp_out = smt.get_empty_model_grid(True).astype(bool)
        temp_wells = wdc_wells.loc[wdc_wells.Zone == z, 'WELL_NO']
        layer, row, col = all_wells.loc[temp_wells, ['layer', 'row', 'col']].values.astype(int).transpose()
        temp_out[layer, row, col] = True
        indexes['{}_well'.format(z.replace(' ', '_'))] = temp_out

    return indexes


def get_cust_indexes():
    """
    get the input indexes for the cust river
    sfr array will be 0 indexed cust reaches in order and -1 where there is no data
    indexes is a dictionary for each sfr_id with the cell flagged to True
    :return: sfr_id_array, indexes  sfr array will be 0 indexed cust reaches in order and -1 where there is no data
    """
    sfr_id_array = smt.shape_file_to_model_array('{}/m_ex_bd_inputs/shp/ordered_cust_reaches.shp'.format(smt.sdp),
                                                 'rid',
                                                 True)
    sfr_id_array[np.isnan(sfr_id_array)] = -1
    sfr_id_array = sfr_id_array.astype(int)
    indexes = {}
    for rid in set(sfr_id_array.flatten()) - {-1}:
        temp = smt.get_empty_model_grid(True).astype(bool)
        temp[0] = np.isclose(sfr_id_array, rid)
        indexes[rid] = temp

    return sfr_id_array, indexes


def get_cust_mapping(base_name, model_ids, recalc=False, recalc_backward_tracking=False):
    """
    set up something to make a dictionary of mappers and return out put, make it save and/or load it so that I don't
    need to keep re-calculating results as netcdf (with geocoding) do both strong/weak forward/backward
    outdata is a dictionary of dictionaries {strong_back: {model_id: packed 3d array}} 3rd dimension is SFR reaches
    losing is a dicionary {model_id: , int8 array) 1 is losing 0 is neither, -1 is gaining
    to pack and unpack bits: see:
    https://stackoverflow.com/questions/5602155/numpy-boolean-array-with-1-bit-entries

    :param base_name: a unique identifier for teh data
    :param model_ids: the model ids to create this for
    :param recalc: bool if True recalc the netcdf
    :param recalc_backward_tracking: bool if True rerun the backward modpath, if flagged true then recalc set to True
    :return outdata, sfr_id_array, unpacked_size, unpacked_shape, losing
    """

    sfr_id_array, indexes = get_cust_indexes()
    outdir = os.path.join(get_base_results_dir('cust', socket.gethostname()), base_name)
    sfr_ids = indexes.keys()
    root_num_part = 4
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = get_base_results_dir('backward', socket.gethostname())
    unpacked_shape = (len(sfr_ids), smt.rows, smt.cols)
    unpacked_size = np.prod(unpacked_shape)

    if recalc_backward_tracking:
        recalc = True
    model_ids = np.atleast_1d(model_ids)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.exists(os.path.join(outdir, 'forward_strong' + '.nc')) and not recalc:
        outdata = {}
        for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
            temp = nc.Dataset(os.path.join(outdir, name + '.nc'))
            temp2 = {e: np.array(temp.variables[e]) for e in
                     set(temp.variables.keys()) - {'latitude', 'longitude', 'sfr_id_array', 'losing', 'model_id'}}
            outdata[name] = temp2
        temp_losing = temp.variables['losing'][:]
        temp_models = temp.variables['model_id'][:]
        losing = {}
        for mid, a in zip(temp_models, temp_losing):
            losing[mid] = a
        return outdata, sfr_id_array, unpacked_size, unpacked_shape, losing

    # forward weak
    print('calculating cust forward weak')
    forward_weaks = []
    f_em_paths = get_forward_emulator_paths(model_ids, True)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        forward_weaks.append(temp)
    forward_weaks = _pack_bits(forward_weaks, model_ids, sfr_ids)

    # forward strong
    print('calculating cust forward strong')
    forward_strongs = []
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        forward_strongs.append(temp)
    forward_strongs = _pack_bits(forward_strongs, model_ids, sfr_ids)

    # backward weak
    print('calculating cust backward weak')
    back_weaks = []
    for mid in model_ids:
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'weak_cust', mid),
                                           mp_name='{}_weak'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=True,
                                           recalc=recalc_backward_tracking)
        back_weaks.append(temp)
    back_weaks = _pack_bits(back_weaks, model_ids, sfr_ids)

    # backward strong
    print('calculating cust backward strong')
    back_strongs = []
    for mid in model_ids:
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'strong_cust', mid),
                                           mp_name='{}_strong'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=False,
                                           recalc=recalc_backward_tracking)
        back_strongs.append(temp)
    back_strongs = _pack_bits(back_strongs, model_ids, sfr_ids)

    losing = {}
    for mid in model_ids:
        temp = smt.get_empty_model_grid().astype(np.int8)
        cbc_path = get_cbc(mid, modflow_dir)
        data = flopy.utils.CellBudgetFile(cbc_path).get_data(kstpkper=(0, 0), text='Stream Leakage', full3D=True)[0][
            0].filled(0)
        limit = 0
        temp[data < -1 * limit] = -1
        temp[data > limit] = 1
        losing[mid] = temp

    # save the data
    print('joining the data')
    outdata = _save_cust_nc(outdir, forward_weaks, forward_strongs, back_strongs, back_weaks, model_ids,
                            root_num_part, sfr_id_array, unpacked_size, unpacked_shape, losing)

    return outdata, sfr_id_array, unpacked_size, unpacked_shape, losing


def create_zones(model_ids, run_name, outdir, root_num_part, indexes, recalc=False, recalc_backward_tracking=False):
    """
    create the zones, load data from pre-run forward models and run and extract the data for backward models
    then amalgamate the data up into useful fashions, also sort out the cust particle tracking problem
    :param model_ids: the model_ids to run this for
    :param outdir: the directory to save the four output netcdfs
    :param root_num_part: the 3 root of the number of particles to put in each backward cell
    :param indexes: a dictionay of ids and boolean zone arrays
    :param recalc: if True recalc the netcdfs
    :param recalc_backward_tracking: if True also re-run backward models (if this is True, recalc will be changed to True
    :return: {forward_weak: open netcdf file}
    """
    if recalc_backward_tracking:
        recalc = True
    model_ids = np.atleast_1d(model_ids)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.exists(os.path.join(outdir, 'forward_strong' + '.nc')) and not recalc:
        outdata = {}
        for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
            outdata[name] = nc.Dataset(os.path.join(outdir, name + '.nc'))
        return outdata

    amalg_weak_forward, amalg_strong_forward, amalg_strong_back, \
    amalg_weak_back, forward_strongs_num_parts = _get_data_for_zones(run_name, model_ids, indexes,
                                                                     root_num_part, recalc_backward_tracking)

    # save the data
    print('joining the data')
    outdata = save_source_nc(outdir=outdir,
                             forward_weak=amalg_weak_forward,
                             forward_strong=amalg_strong_forward,
                             backward_strong=amalg_strong_back,
                             backward_weak=amalg_weak_back,
                             model_ids=model_ids,
                             root_num_part=root_num_part,
                             forward_part_nums=forward_strongs_num_parts)

    return outdata


def save_source_nc(outdir, forward_weak, forward_strong, backward_strong, backward_weak, model_ids, root_num_part,
                   forward_part_nums):
    """
    save the thing as a netcdf file
    dimensions for data(direction, source behavior, amalg_type, row, col) variables of location ids
    :param outdir: directory to save the four files in
    :param forward_weak: the data for the forward weak runs
    :param forward_strong: the data for the forward strong runs
    :param backward_strong: the data for the backward strong runs
    :param backward_weak: the data for the backward weak runs
    :param model_ids: a list of model ids
    :param root_num_part: the 3 root of the number of particles released in each backward cell
    :param forward_part_nums: not used now, but the list of the forward particle runs (here incase I want to save it later)
    :return: {backward_strong: netcdf file}
    """
    # some assertions
    assert all([isinstance(e, dict) for e in
                [forward_strong, forward_weak, backward_strong, backward_weak]]), 'all data must be dict'
    temp = set(forward_weak.keys()) == set(forward_strong.keys()) == set(backward_weak.keys()) == set(
        backward_strong.keys())
    assert temp, 'all inputs must have the same keys'
    no_flow = smt.get_no_flow(0)
    amalg_types = list(set([e.split('~')[-1] for e in forward_weak.keys()]))
    sites = list(set([e.split('~')[0] for e in forward_weak.keys()]))

    outdata = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        outfile = nc.Dataset(os.path.join(outdir, name + '.nc'), 'w')
        x, y = smt.get_model_x_y(False)
        # create dimensions
        outfile.createDimension('latitude', len(y))
        outfile.createDimension('longitude', len(x))

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

        # location add the data
        for site in eval(name).keys():
            temp_var = outfile.createVariable(site, 'u1',  # this is u1(np.uint8) because we are only running 165 sims?
                                              ('latitude', 'longitude'),
                                              fill_value=0)
            temp_var.setncatts({'units': 'bool or number of realisations',
                                'long_name': site,
                                'missing_value': 0,
                                'comments': 'number of realisations that show particles from a given cell'})

            t = eval(name)[site].astype(np.uint8)
            t[~np.isclose(no_flow, 1)] = 0
            temp_var[:] = t

        outfile.description = ('source zones for single sources')
        outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.backward_num_parts = root_num_part ** 3
        outfile.model_ids = ' '.join(model_ids)
        outfile.close()

        outdata[name] = nc.Dataset(os.path.join(outdir, name + '.nc'))  # reload to avoid write behavior
    return outdata


def _save_cust_nc(outdir, forward_weak, forward_strong, backward_strong, backward_weak, model_ids, root_num_part,
                  sfr_ids_array, unpacked_size, unpacked_shape, losing):
    """
    save the cust id  as a netcdf file
    dimensions for data(direction, source behavior, amalg_type, row, col) variables of location ids

    :param outdir: directory to save the four files in
    :param forward_weak: the data for the forward weak runs
    :param forward_strong: the data for the forward strong runs
    :param backward_strong: the data for the backward strong runs
    :param backward_weak: the data for the backward weak runs
    :param model_ids: a list of model ids
    :param root_num_part: the 3 root of the number of particles released in each backward cell
    :param sfr_ids_array: an array identifing the sfr ids
    :param unpacked_size: the unpacked size of the packed 3d array
    :param unpacked_shape:  the unpacked shape of the packed 3d array (sfr_ids, rows, cols)
    :param losing: {model_id: boolean 2d array for gaining and/or losing reach} True is losing
    :return: the packed data {strong_back: {model_id: packed 3d array}}
    """
    # some assertions
    assert all([isinstance(e, dict) for e in
                [forward_strong, forward_weak, backward_strong, backward_weak]]), 'all data must be dict'
    temp = set(forward_weak.keys()) == set(forward_strong.keys()) == set(backward_weak.keys()) == set(
        backward_strong.keys())
    assert temp, 'all inputs must have the same keys'

    outdata = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        outfile = nc.Dataset(os.path.join(outdir, name + '.nc'), 'w')
        x, y = smt.get_model_x_y(False)
        # create dimensions
        outfile.createDimension('latitude', len(y))
        outfile.createDimension('longitude', len(x))
        outfile.createDimension('packed_dim')
        outfile.createDimension('model_id')

        # create variables

        mid = outfile.createVariable('model_id', str, ('model_id',))
        mid[:] = model_ids

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

        cust_array = outfile.createVariable('sfr_id_array', int, ('latitude', 'longitude'), fill_value=-1)
        cust_array.setncatts({'units': 'None',
                              'long_name': 'indexes for sfr array',
                              'missing_value': -1,
                              })
        cust_array[:] = sfr_ids_array

        losing_array = np.concatenate([losing[e][np.newaxis] for e in model_ids])
        los = outfile.createVariable('losing', int, ('model_id', 'latitude', 'longitude'), fill_value=0)
        los.setncatts({'units': 'None',
                       'long_name': 'losing/gaining reaches -1 gaining, 1 is losing, 0 is neither/nodata',
                       'missing_value': 0,
                       })
        los[:] = losing_array

        # location add the data
        for mid in eval(name).keys():
            temp_var = outfile.createVariable(mid, 'u1',
                                              ('packed_dim',),
                                              fill_value=0)
            temp_var.setncatts({'units': 'packed boolean',
                                'long_name': mid,
                                'missing_value': 0,
                                'unpacked_size': unpacked_size,
                                'unpacked_shape': str(unpacked_shape)})

            t = eval(name)[mid]
            temp_var[:] = t

        outfile.description = ('source zones for the cust')
        outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.backward_num_parts = root_num_part ** 3
        outfile.model_ids = ' '.join(model_ids)
        outfile.upacked_size = unpacked_size
        outfile.upacked_shape = str(unpacked_shape)  # the tuple can be called by using eval
        outfile.close()

        temp = nc.Dataset(os.path.join(outdir, name + '.nc'))
        temp2 = {e: np.array(temp.variables[e]) for e in
                 set(temp.variables.keys()) - {'latitude', 'longitude', 'sfr_id_array', 'losing', 'model_id'}}
        outdata[name] = temp2
    return outdata


def _pack_bits(data, model_ids, sfr_ids):
    """
    packs and reformats the data
    :param data: the data [{sfr_id[0]:array}, ..., len(model_ids[-1])] _see forward_weaks as example
    :param model_ids: a list of model ids
    :param sfr_ids: a list of sfr_ids
    :return: outdata {model_id: packed 3d array}
    """
    # change the outputs of the define sources to boolean and pack into bits
    assert isinstance(data, list), 'data must be dict'
    assert len(data) == len(model_ids), 'len data must equal len of model ids'
    assert (all([isinstance(e, dict) for e in data])), 'all entries must be dict'
    assert (all([np.in1d(sfr_ids, e.keys()).all() for e in data])), 'all sfr_ids must be present in data'

    outdata = {}
    for i, mid in enumerate(model_ids):
        temp = np.concatenate([data[i][e][np.newaxis] for e in sfr_ids], 0).astype(bool)
        outdata[mid] = np.packbits(temp)
    return outdata


def _get_data_for_zones(run_name, model_ids, indexes, root_num_part, recalc_backward_tracking):
    """
    get and amalgamate up the data for the zone delination
    :param run_name: the name to call this run of backward models (to prevent overwrite)
    :param model_ids: list of model ids
    :param indexes: dictionary of boolean arrays for the targets
    :param root_num_part: the cubic root of the number of particles to release in each backward modpath cell
    :param recalc_backward_tracking: bool, if true re-run the backward partical tracking
    :return: amalg_weak_forward, amalg_strong_forward, amalg_strong_back, amalg_weak_back, forward_strongs_num_parts
    """
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = os.path.join(get_base_results_dir('backward', socket.gethostname()), run_name)

    cust_data = get_cust_mapping(run_name, model_ids)

    # backward weak
    print('calculating backward weak\n\n')
    back_weaks = []
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i+1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'weak', mid),
                                           mp_name='{}_weak'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=True,
                                           recalc=recalc_backward_tracking)
        back_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_weak_back = _amalg_backward(back_weaks, indexes, cust_data, model_ids, 'backward_weak')

    # backward strong
    print('calculating backward strong\n\n')
    back_strongs = []
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i+1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'strong', mid),
                                           mp_name='{}_strong'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=False,
                                           recalc=recalc_backward_tracking)
        back_strongs.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_strong_back = _amalg_backward(back_strongs, indexes, cust_data, model_ids, 'backward_strong')

    # forward weak
    print('calculating forward weak\n\n')
    forward_weaks = []
    forward_weaks_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, True)

    for i, path in enumerate(f_em_paths.values()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i+1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        temp2 = np.loadtxt(path[1].replace('_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_weaks_num_parts.append(temp2)
        forward_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_weak_forward = _amalg_forward(forward_weaks, indexes, cust_data, model_ids, 'forward_weak')

    # forward strong
    print('calculating forward strong\n\n')
    forward_strongs = []
    forward_strongs_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)

    for i, path in enumerate(f_em_paths.values()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i+1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        forward_strongs.append(temp)
        temp2 = np.loadtxt(path[1].replace('_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_strongs_num_parts.append(temp2)
    # amalgamate data from realsiations
    amalg_strong_forward = _amalg_forward(forward_strongs, indexes, cust_data, model_ids, 'forward_strong')

    return amalg_weak_forward, amalg_strong_forward, amalg_strong_back, amalg_weak_back, forward_strongs_num_parts


def _amalg_forward(data, indexes, sfr_data, model_ids, run_name):
    """
    :param data: the list of dictionaries for the data
    :param indexes: the indexes going forward
    :param model_ids: the list of model ids in teh same order as data
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return: {id: map} id should match amalg backward
    """
    out = {}
    for name in indexes.keys():
        temp = _get_all_avalible_data(data, name)
        _add_data_variations(out, temp, name, sfr_data, model_ids, run_name)

    return out


def _amalg_backward(data, indexes, sfr_data, model_ids, run_name):
    """

    :param data: the list of dictionaries for the data
    :param indexes: the indexes going forward
    :param model_ids: the list of model ids in teh same order as data
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return: {id: map} id should match amalg forward
    """
    out = {}
    for name in indexes.keys():
        temp = _get_all_avalible_data(data, name)
        _add_data_variations(out, temp, name, sfr_data, model_ids, run_name)

    return out

def _get_all_avalible_data(data,key):
    temp = []
    for e in data:
        try:
            temp.append(e[key])
        except:
            pass
    return temp


def _add_data_variations(out, org_arrays, name, sfr_data, model_ids, run_name):
    """
    add amagimated data
    :param out: dictionay to add the amalgamated data to
    :param org_arrays: a list of the original arrays
    :param name: base name to call the variable
    :param sfr_data: the outputs (as a tuple) for get cust mapping
    :param model_ids: the list of model ids in teh same order as org arrays
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return:
    """

    sfr_data, sfr_id_array, unpacked_size, unpacked_shape, losing = sfr_data
    bool_array = [e > 0 for e in org_arrays]
    out['{}_all'.format(name)] = np.all(bool_array, axis=0).astype(np.uint8)
    out['{}_number'.format(name)] = np.sum(bool_array, axis=0).astype(np.uint8)

    # add upstream cust influance
    bool_array_wcust = []
    for mid, temp_bool_array in zip(model_ids, bool_array):
        temp_sfr_id_array = deepcopy(sfr_id_array)
        temp_sfr_id_array[losing[mid] < 0] = -1

        if not (temp_sfr_id_array[temp_bool_array] >= 0).any():
            bool_array_wcust.append(temp_bool_array)
            continue

        temp_ids = np.array(list(set(temp_sfr_id_array[temp_bool_array]) - {-1})).astype(int)
        # get and unpack the array
        temp = np.unpackbits(sfr_data[run_name][mid])[:unpacked_size].reshape(unpacked_shape).astype(bool)
        temp = temp[:temp_ids.max() + 1].sum(axis=0)
        temp += temp_bool_array

        bool_array_wcust.append(temp.astype(bool))

    out['{}_all_cust'.format(name)] = np.all(bool_array_wcust, axis=0).astype(np.uint8)
    out['{}_number_cust'.format(name)] = np.sum(bool_array_wcust, axis=0).astype(np.uint8)


def run_single_source_zones(recalc=False, recalc_backward_tracking=False):
    """
    run the singel source zones... final wrapper
    :param recalc: bool if True recalculate the netcdf file
    :param recalc_backward_tracking: bool if True recalculate the modpath particle tracking simulations
    :return:
    """
    indexes = create_single_zone_indexs()
    base_outdir = r"C:\mh_waimak_models\single_source_zones"
    print('running for AshOpt')
    outdir = os.path.join(base_outdir, 'AshOpt')
    create_zones(model_ids=['AshOpt'], run_name='AshOpt_single_sources',
                 outdir=outdir, root_num_part=3,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(outdir)

    print('running for 165 models')
    outdir = os.path.join(base_outdir, 'stocastic set')
    stocastic_model_ids = get_stocastic_set()
    create_zones(model_ids=stocastic_model_ids,
                 run_name='stocastic_set_single_sources',
                 outdir=outdir, root_num_part=3,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(outdir)


def split_netcdfs(indir):
    """
    split up the netcdfs into netcdfs for the sites
    :param indir: the dir with the 3 netcdfs  the new files are put in a file labeled individuals
    :return:
    """
    print('splitting netcdf')
    outdir = os.path.join(indir, 'individual_netcdfs')
    if not os.path.exists(os.path.join(indir, 'individual_netcdfs')):
        os.makedirs(outdir)

    data = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        data[name] = nc.Dataset(os.path.join(indir, name + '.nc'))

    # get list of variables (assmue all are teh same) and list of base variables
    variables = list(set(data['forward_strong'].variables.keys()) - {'crs', 'latitude', 'longitude'})
    base_variables = list(set([e.replace('_number', '').replace('_all', '').replace('_cust', '') for e in variables]))

    for bv in base_variables:
        outfile = nc.Dataset(os.path.join(outdir, bv + '.nc'), 'w')
        x, y = smt.get_model_x_y(False)
        # create dimensions
        outfile.createDimension('latitude', len(y))
        outfile.createDimension('longitude', len(x))

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

        # location add the data
        for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
            for suffix in ['number', 'number_cust']:
                temp_var = outfile.createVariable('{}_{}_{}'.format(name[0:4], name.split('_')[-1][0:2], suffix), float,
                                                  ('latitude', 'longitude'),
                                                  fill_value=np.nan)
                temp_var.setncatts({'units': 'bool or number of realisations',
                                    'long_name': '{}_{}_{}'.format(name, bv, suffix),
                                    'missing_value': np.nan,
                                    'comments': 'number of particles from a given cell'})

                t = data[name].variables['{}_{}'.format(bv, suffix)][:]
                t = t.filled(0).astype(float)
                t[np.isclose(t, 0)] = np.nan
                temp_var[:] = t

        outfile.description = ('source zones for single sources')
        outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.close()


if __name__ == '__main__':
    run_single_source_zones()