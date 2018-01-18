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
import traceback
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.source_delmination.source_raster import \
    define_source_from_backward, define_source_from_forward, get_modeflow_dir_for_source, get_base_results_dir, \
    get_cbc, get_forward_emulator_paths


# set up all single zone delination maps (e.g. streams/ water supply wells)
# todo how to handle the cust/+- ashley problem? manually? that would be a pain in the ass I;m not sure what to do here...
# I could run one for all of the cust and add the particles... would mess up my numbering normalisation...

def create_single_zone_indexs():  # todo
    """
    create an index for all zones
    :return: {id: index}
    """

    # todo the below is just a tester
    index = smt.get_empty_model_grid(True).astype(bool)
    index = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\rough_chch.shp".format(smt.sdp), 'Id', True)[
        np.newaxis]
    index2 = np.isfinite(index).repeat(11, axis=0)
    index2[1:, :, :] = False
    index1 = np.full((smt.layers, smt.rows, smt.cols), False)
    index1[6] = np.isfinite(index[0])
    indexes = {'layer0_chch': index2, 'layer7_chch': index1}
    return indexes


def create_zones(model_ids, outdir, root_num_part, recalc=False, recalc_backward_tracking=False):
    """
    set up something to make a dictionary of mappers and return out put, make it save and/or load it so that I don't
    need to keep re-calculating results as netcdf (with geocoding) do both strong/weak forward/backward
    :return:
    """
    if recalc_backward_tracking:
        recalc = True
    model_ids = np.atleast_1d(model_ids)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if os.path.exists(os.path.join(outdir, 'forward_strong' + '.nc')) and not recalc:
        outdata = {}
        for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
            outdata[name] = nc.Dataset(os.path.join(outdir,name+'.nc'))
        return outdata

    indexes = create_single_zone_indexs()
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = get_base_results_dir('backward', socket.gethostname())


    # forward weak
    print('calculating forward weak')
    forward_weaks = []
    forward_weaks_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, True)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        temp2 = np.loadtxt(path[1].replace('_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_weaks_num_parts.append(temp2)
        forward_weaks.append(temp)
    # amalgamate data from realsiations mean and sd?
    amalg_weak_forward = _amalg_forward(forward_weaks, indexes)

    # forward strong
    print('calculating forward strong')
    forward_strongs = []
    forward_strongs_num_parts = []
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)
    for path in f_em_paths.values():
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        forward_strongs.append(temp)
        temp2 = np.loadtxt(path[1].replace('_bnd_type.txt', '_num_parts.txt'))  # load number of particles
        forward_strongs_num_parts.append(temp2)
    # amalgamate data from realsiations
    amalg_strong_forward = _amalg_forward(forward_strongs, indexes)

    # backward weak
    print('calculating backward weak')
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
    amalg_weak_back = _amalg_backward(back_weaks, indexes)

    # backward strong
    print('calculating backward strong')
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
    amalg_strong_back = _amalg_backward(back_strongs, indexes)

    # save the data
    print('joining the data')
    outdata = _save_source_nc(outdir=outdir,
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


def _save_source_nc(outdir, forward_weak, forward_strong, backward_strong, backward_weak, model_ids, root_num_part,
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
    no_flow = smt.get_no_flow(0)
    amalg_types = list(set([e.split('~')[-1] for e in forward_weak.keys()]))
    sites = list(set([e.split('~')[0] for e in forward_weak.keys()]))

    outdata = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        outfile = nc.Dataset(os.path.join(outdir, name+'.nc'), 'w')
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
            temp_var = outfile.createVariable(site, float,
                                              ('latitude', 'longitude'),
                                              fill_value=np.nan)
            temp_var.setncatts({'units': 'fraction',
                                'long_name': site,
                                'missing_value': np.nan,
                                'comments': 'number of particles from a given cell'})

            t = eval(name)[site].astype(float)
            t[~np.isclose(no_flow, 1)] = np.nan
            temp_var[:] = t

        outfile.description = ('source zones for single sources')
        outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.backward_num_parts = root_num_part ** 3
        outfile.model_ids = ' '.join(model_ids)
        outfile.close()

        outdata[name] = nc.Dataset(os.path.join(outdir, name + '.nc'))  # reload to avoid write behavior
    return outdata


def _amalg_forward(data, indexes):
    out = {}
    for name in indexes.keys():
        temp = [e[name] for e in data]
        _add_data_variations(out, temp, name)

    return out


def _amalg_backward(data, indexes):
    """

    :param data: the list of dictionaries for the data
    :param root_num_part: the root number of particles generated
    :param indexes: the indexes going forward
    :return: {id: map} id should match amalg forward
    """
    out = {}
    for name in indexes.keys():
        temp = [e[name] for e in data]
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
    if False: # I did not find these terribly useful at present
        out['{}_u'.format(name)] = np.mean(temp, axis=0)
        out['{}_sd'.format(name)] = np.std(temp, axis=0)
        out['{}_sum'.format(name)] = np.sum(temp, axis=0)
        out['{}_min'.format(name)] = np.min(temp, axis=0)
        out['{}_max'.format(name)] = np.min(temp, axis=0)
        out['{}_5th'.format(name)] = np.percentile(temp, 5, axis=0)
        out['{}_25th'.format(name)] = np.percentile(temp, 25, axis=0)
        out['{}_50th'.format(name)] = np.percentile(temp, 50, axis=0)
        out['{}_75th'.format(name)] = np.percentile(temp, 75, axis=0)
        out['{}_95th'.format(name)] = np.percentile(temp, 95, axis=0)
    temp2 = [e > 0 for e in temp]
    out['{}_all'.format(name)] = np.all(temp2, axis=0).astype(float)
    out['{}_any'.format(name)] = np.any(temp2, axis=0).astype(float)
    out['{}_number'.format(name)] = np.sum(temp2, axis=0).astype(float)


if __name__ == '__main__':
    create_zones(model_ids=['NsmcBase', 'AshOpt'], outdir=r"C:\Users\matth\Downloads\test_zone_delin",
                 root_num_part=3, recalc=True, recalc_backward_tracking=False)
