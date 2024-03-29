# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/01/2018 10:26 AM
"""

from __future__ import division
from single_zone_delination import get_cust_mapping
import numpy as np
import netCDF4 as nc
import os
import socket
import sys
import datetime
import pandas as pd
from copy import deepcopy
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.source_raster import \
    define_source_from_backward, define_source_from_forward, get_modeflow_dir_for_source, get_base_results_dir, \
    get_forward_emulator_paths
from waimak_extended_boundry.model_run_tools.model_bc_data.all_well_layer_col_row import get_all_well_row_col
from waimak_extended_boundry.model_run_tools import get_stocastic_set
import gc


def create_private_wells_indexes(num_iterations=1, iteration=0):
    """
    get the private well indexes
    :param num_iterations: the number of itterations to break this up into
    :param iteration: which iteration
    :return: out_private_wells, indexes
    """
    all_wells = get_all_well_row_col()

    private_wells = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\PrivateWellZones.csv",
        index_col=0)
    private_wells = pd.merge(private_wells, all_wells.loc[:, ['layer', 'row', 'col']], right_index=True,
                             left_index=True)
    private_wells = private_wells.dropna()
    private_wells = pd.merge(private_wells, pd.DataFrame(all_wells.loc[:, 'depth']), how='left', left_index=True,
                             right_index=True)
    idx = private_wells.depth > 50
    private_wells.loc[idx, 'zone_2'] = private_wells.loc[idx, 'Zone_1'] + '_deep'
    private_wells.loc[~idx, 'zone_2'] = private_wells.loc[~idx, 'Zone_1'] + '_shallow'

    chunk_size = len(private_wells) // num_iterations + 1
    if num_iterations == iteration + 1:  # if it is the last iteration then use +1 to grab the last index
        extra = 1
    else:
        extra = 0
    start = iteration * chunk_size
    end = start + chunk_size + extra
    out_private_wells = private_wells.iloc[start:end]

    indexes = {}
    for name, k, i, j in out_private_wells.loc[:, ['layer', 'row', 'col']].astype(int).itertuples(True, None):
        temp = smt.get_empty_model_grid(True).astype(bool)
        temp[k, i, j] = True
        indexes[name] = temp

    return out_private_wells, indexes


def create_amalgimated_source_protection_zones(model_ids, run_name, outdir, num_it, recalc=False,
                                               recalc_backward_tracking=False):
    """
    calculates teh source area for all the individaul wells (or loads it) and then aggrigates the data across the
    different zones new amalgimated netcdfs are placed in a amalgimate folder in outdir  all arguments are passed
    directly to create_zones
    :param model_ids:
    :param run_name:
    :param outdir:
    :param recalc:
    :param recalc_backward_tracking:
    :return:
    """
    private_wells, indexes = create_private_wells_indexes()
    root_num_part = 3
    single_site_data = create_zones(model_ids, run_name, outdir, root_num_part, num_it,
                                    recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)

    # the amalgimation must begin
    for key in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        data = single_site_data[key]
        outdata = {}
        sites = {}
        for site in set(private_wells.Zone_1):
            well_nums = list(private_wells.loc[private_wells.Zone_1 == site].index)
            sites[site] = well_nums
        for site in set(private_wells.zone_2):
            well_nums = list(private_wells.loc[private_wells.zone_2 == site].index)
            sites[site] = well_nums
        sites['all_wells'] = list(private_wells.index)
        sites['all_wells_shallow'] = list(private_wells.loc[private_wells.depth <= 50].index)
        sites['all_wells_deep'] = list(private_wells.loc[private_wells.depth > 50].index)

        for site, well_nums in sites.items():
            any_array = smt.get_empty_model_grid().astype(int)
            all_array = smt.get_empty_model_grid().astype(int)
            number_array = smt.get_empty_model_grid().astype(int)
            any_array_cust = smt.get_empty_model_grid().astype(int)
            all_array_cust = smt.get_empty_model_grid().astype(int)
            number_array_cust = smt.get_empty_model_grid().astype(int)
            for well in well_nums:
                well = well.replace('/', '_')
                temp_all = np.array(data.variables['{}_all'.format(well)])
                temp_all_cust = np.array(data.variables['{}_all_cust'.format(well)])
                temp_number = np.array(data.variables['{}_number'.format(well)])
                temp_number_cust = np.array(data.variables['{}_number_cust'.format(well)])

                number_array += temp_number.astype(int)
                number_array_cust += temp_number_cust.astype(int)

                any_array += (temp_number > 0).astype(int)
                any_array_cust += (temp_number_cust > 0).astype(int)

                all_array += temp_all.astype(int)
                all_array_cust += temp_all_cust.astype(int)
            outdata['{}_all'.format(site)] = all_array
            outdata['{}_all_cust'.format(site)] = all_array_cust
            outdata['{}_any'.format(site)] = any_array
            outdata['{}_any_cust'.format(site)] = any_array_cust
            outdata['{}_number'.format(site)] = number_array
            outdata['{}_number_cust'.format(site)] = number_array_cust

        save_source_nc(outdir, 'amalgimated_{}'.format(key), outdata, model_ids, root_num_part, True, dtype=np.uint64)


def create_zones(model_ids, run_name, outdir, root_num_part, num_iterations=1, recalc=False,
                 recalc_backward_tracking=False):
    """
    create the zones, load data from pre-run forward models and run and extract the data for backward models
    then amalgamate the data up into useful fashions, also sort out the cust particle tracking problem
    :param model_ids: the model_ids to run this for
    :param outdir: the directory to save the four output netcdfs
    :param root_num_part: the 3 root of the number of particles to put in each backward cell
    :param num_iterations: int the number of iterations to run this thing in
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

    modflow_dir = get_modeflow_dir_for_source()

    cust_data = get_cust_mapping(run_name, model_ids)

    for it in range(num_iterations):
        print('\n\n##################################################')
        print('starting well set {} of {}'.format(it + 1, num_iterations))
        print('###################################################\n\n')
        private_wells, indexes = create_private_wells_indexes(num_iterations, it)
        backward_dir = os.path.join(get_base_results_dir('backward', socket.gethostname()),
                                    '{}_it{}'.format(run_name, it))
        if it == 0:
            first = True
        else:
            first = False

        # forward weak
        print('calculating forward weak\n\n')
        forward_weaks = []
        f_em_paths = get_forward_emulator_paths(model_ids, True)
        for i, path in enumerate(f_em_paths.values()):
            print('{}, {} of {}'.format(os.path.basename(path[0]), i + 1, len(f_em_paths)))
            temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes,
                                              return_packed_bits=True)
            forward_weaks.append(temp)

        # amalgamate data from realsiations mean and sd?
        amalg_weak_forward = _amalg_forward(forward_weaks, indexes, cust_data, model_ids, 'forward_weak')
        outpath = save_source_nc(outdir=outdir,
                                 name='forward_weak',
                                 data=amalg_weak_forward,
                                 model_ids=model_ids,
                                 root_num_part=root_num_part,
                                 new_file=first)
        gc.collect()

        # forward strong
        print('calculating forward strong\n\n')
        forward_strongs = []
        f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)
        for i, path in enumerate(f_em_paths.values()):
            print('{}, {} of {}'.format(os.path.basename(path[0]), i + 1, len(f_em_paths)))
            temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes,
                                              return_packed_bits=True)
            forward_strongs.append(temp)
        # amalgamate data from realsiations
        amalg_strong_forward = _amalg_forward(forward_strongs, indexes, cust_data, model_ids, 'forward_strong')
        outpath = save_source_nc(outdir=outdir,
                                 name='forward_strong',
                                 data=amalg_strong_forward,
                                 model_ids=model_ids,
                                 root_num_part=root_num_part,
                                 new_file=first)
        gc.collect()

        # backward weak
        print('calculating backward weak\n\n')
        back_weaks = []
        for i, mid in enumerate(model_ids):
            print('model: {}, {} of {}'.format(mid, i + 1, len(model_ids)))
            temp = define_source_from_backward(indexes,
                                               mp_ws=os.path.join(backward_dir, 'weak', mid),
                                               mp_name='{}_weak'.format(mid),
                                               cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                               root3_num_part=root_num_part, capt_weak_s=True,
                                               recalc=recalc_backward_tracking,
                                               return_packed_bits=True)
            back_weaks.append(temp)
        # amalgamate data from realsiations mean and sd?
        # there are some wells are missing, my guess is they got stranded (which accounts for most (e.g. dry wells)
        # and/or the exit via a non top layer flux
        amalg_weak_back = _amalg_backward(back_weaks, indexes, cust_data, model_ids,
                                          'backward_weak')  # there are some wells are missing, my guess is they got stranded (which accounts for most (e.g. dry wells) and/or the exit via a non top layer flux
        outpath = save_source_nc(outdir=outdir,
                                 name='backward_weak',
                                 data=amalg_weak_back,
                                 model_ids=model_ids,
                                 root_num_part=root_num_part,
                                 new_file=first)

        gc.collect()

        # backward strong
        print('calculating backward strong\n\n')
        back_strongs = []
        for i, mid in enumerate(model_ids):
            print('model: {}, {} of {}'.format(mid, i + 1, len(model_ids)))
            temp = define_source_from_backward(indexes,
                                               mp_ws=os.path.join(backward_dir, 'strong', mid),
                                               mp_name='{}_strong'.format(mid),
                                               cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                               root3_num_part=root_num_part, capt_weak_s=False,
                                               recalc=recalc_backward_tracking, return_packed_bits=True)
            back_strongs.append(temp)
        # amalgamate data from realsiations mean and sd?
        amalg_strong_back = _amalg_backward(back_strongs, indexes, cust_data, model_ids, 'backward_strong')
        outpath = save_source_nc(outdir=outdir,
                                 name='backward_strong',
                                 data=amalg_strong_back,
                                 model_ids=model_ids,
                                 root_num_part=root_num_part,
                                 new_file=first)

        gc.collect()

    outdata = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        outdata[name] = nc.Dataset(os.path.join(outdir, name + '.nc'))
    return outdata


def save_source_nc(outdir, name, data, model_ids, root_num_part, new_file=False, dtype=np.uint8):
    """
    save the thing as a netcdf file
    dimensions for data(direction, source behavior, amalg_type, row, col) variables of location ids
    :param outdir: directory to save the four files in
    :param name: one of the standard names (ends up being the file name
    :param data: the data for the source area
    :param model_ids: a list of model ids
    :param root_num_part: the 3 root of the number of particles released in each backward cell
    :return: path to netcdf file
    """
    # some assertions
    assert isinstance(data, dict), 'data must be dict'
    no_flow = smt.get_no_flow(0)

    outpath = os.path.join(outdir, name + '.nc')
    if new_file:
        outfile = nc.Dataset(outpath, 'w')
        outfile.set_fill_off()

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

        lat = outfile.createVariable('latitude', 'f8', ('latitude',))
        lat.setncatts({'units': 'NZTM',
                       'long_name': 'latitude',
                       'standard_name': 'projection_y_coordinate'})
        lat[:] = y

        lon = outfile.createVariable('longitude', 'f8', ('longitude',))
        lon.setncatts({'units': 'NZTM',
                       'long_name': 'longitude',
                       'standard_name': 'projection_x_coordinate'})
        lon[:] = x
    else:
        outfile = nc.Dataset(outpath, 'a')

    # location add the data

    for site in data.keys():
        temp_var = outfile.createVariable(site.replace('/', '_'), dtype, ('latitude', 'longitude'))
        temp_var.setncatts({'units': 'bool or number of realisations',
                            'long_name': site,
                            'comments': 'number of particles from a given cell'})

        t = data[site].astype(dtype)
        t[~np.isclose(no_flow, 1)] = 0
        temp_var[:] = t

    outfile.description = ('source zones for single sources')
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'script: {}'.format(sys.argv[0])
    outfile.backward_num_parts = root_num_part ** 3
    outfile.model_ids = ' '.join(model_ids)
    outfile.close()

    return outpath


def _amalg_forward(data, indexes, sfr_data, model_ids, run_name):
    """
    :param data: the list of dictionaries for the data
    :param root_num_part: the root number of particles generated
    :param indexes: the indexes going forward
    :param model_ids: the list of model ids in teh same order as data
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return: {id: map} id should match amalg backward
    """
    out = {}
    for name in indexes.keys():
        temp = []
        for e in data:  # set all False arrays for those wells who never intersect particles
            try:
                temp.append(e[name])
            except KeyError:
                temp.append(np.zeros((16608), dtype=np.uint8))

        _add_data_variations(out, temp, name, sfr_data, model_ids, run_name)

    return out


def _amalg_backward(data, indexes, sfr_data, model_ids, run_name):
    """

    :param data: the list of dictionaries for the data
    :param root_num_part: the root number of particles generated
    :param indexes: the indexes going forward
    :param model_ids: the list of model ids in teh same order as data
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return: {id: map} id should match amalg forward
    """
    out = {}
    for name in indexes.keys():
        temp = []
        for e in data:  # set all False arrays for those wells who's particles get stranded
            try:
                temp.append(e[name])
            except KeyError:
                temp.append(np.zeros((16608), dtype=np.uint8))

        _add_data_variations(out, temp, name, sfr_data, model_ids, run_name)

    return out


def _add_data_variations(out, org_arrays_packed, name, sfr_data, model_ids, run_name):
    """
    add amagimated data
    :param out: dictionay to add the amalgamated data to
    :param org_arrays: a list of the original arrays (packed_bits)
    :param name: base name to call the variable
    :param sfr_data: the outputs (as a tuple) for get cust mapping
    :param model_ids: the list of model ids in teh same order as org arrays
    :param run_name: one of ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
    :return:
    """
    bool_array = [np.unpackbits(e)[:132860].reshape((smt.rows, smt.cols)).astype(bool) for e in org_arrays_packed]
    sfr_data, sfr_id_array, unpacked_size, unpacked_shape, losing = sfr_data

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


def run_multiple_source_zones(num_it_stocastic, recalc=False, recalc_backward_tracking=False):
    """
    the raper to
    :param num_it_stocastic: the number of iterations for the stocastic
    :param recalc: bool if true recalculate the netcdf
    :param recalc_backward_tracking: bool, if true re-run the modpath simulations
    :return:
    """
    base_outdir = r"D:\mh_waimak_models\private_domestic_supply"
    print('running for AshOpt')
    create_amalgimated_source_protection_zones(model_ids=['AshOpt'], run_name='AshOpt_private_wells',
                                               outdir=os.path.join(base_outdir, 'AshOpt'), num_it=1,
                                               recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(os.path.join(base_outdir, 'AshOpt'))

    print('running for 165 models')
    stocastic_model_ids = get_stocastic_set()
    create_amalgimated_source_protection_zones(model_ids=stocastic_model_ids,
                                               run_name='stocastic_set_private_wells',
                                               outdir=os.path.join(base_outdir, 'stocastic set'),
                                               num_it=num_it_stocastic,
                                               recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(os.path.join(base_outdir, 'stocastic set'))


def split_netcdfs(indir):
    """
    split the output netcdfs into netcdfs by variables
    :param indir: the dir with the 3 netcdfs  the new files are put in a file labeled individuals
    :return:
    """
    print('splitting netcdf')
    outdir = os.path.join(indir, 'individual_netcdfs')
    if not os.path.exists(os.path.join(indir, 'individual_netcdfs')):
        os.makedirs(outdir)

    data = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        data[name] = nc.Dataset(os.path.join(indir, 'amalgimated_' + name + '.nc'))

    # get list of variables (assmue all are teh same) and list of base variables
    variables = list(set(data['forward_strong'].variables.keys()) - {'crs', 'latitude', 'longitude'})
    base_variables = list(
        set([e.replace('_number', '').replace('_all', '').replace('_cust', '').replace('_any', '') for e in variables]))

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
            for suffix in ['number', 'number_cust', 'all', 'all_cust', 'any', 'any_cust']:
                temp_var = outfile.createVariable('{}_{}_{}'.format(name[0:4], name.split('_')[-1][0:2], suffix), float,
                                                  ('latitude', 'longitude'),
                                                  fill_value=np.nan)
                temp_var.setncatts({'units': 'bool or number of realisations',
                                    'long_name': '{}_{}_{}'.format(name, bv, suffix),
                                    'missing_value': np.nan,
                                    'comments': 'number of particles from a given cell'})

                t = data[name].variables['{}_{}'.format(bv, suffix)][:]
                t = t.astype(float)
                t[np.isclose(t, 0)] = np.nan
                temp_var[:] = t

        outfile.description = ('source zones for aggregated sources. '
                               'any: sum of wells where any model showed the cell in a source zone '
                               'all: sum of wells where all models showed the cell in a source zone '
                               'number: num of the number of models which showed the cell in a source zone '
                               '(also summed over wells) '
                               '_cust: denotes particle tracking through the cust river')
        outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.close()


if __name__ == '__main__':
    run_multiple_source_zones(3, recalc=False, recalc_backward_tracking=False)
