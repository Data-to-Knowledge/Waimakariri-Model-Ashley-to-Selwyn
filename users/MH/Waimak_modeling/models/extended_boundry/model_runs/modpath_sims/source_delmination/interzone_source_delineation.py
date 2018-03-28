"""
Author: matth
Date Created: 27/03/2018 7:35 PM
"""

from __future__ import division
from core import env
import os
import numpy as np
import netCDF4 as nc
import socket
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set
from single_zone_delination import get_modeflow_dir_for_source, get_base_results_dir, get_cust_mapping, \
    define_source_from_backward, get_cbc, get_forward_emulator_paths, define_source_from_forward
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from datetime import datetime
import sys
from copy import deepcopy
import geopandas as gpd
import itertools


base_receptors_path = env.sci(
    r"Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\shp\interzone_receptors.shp")

def make_shapefiles(outdir):
    # make shapefiles for viewing
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    base_sites = get_base_sites()
    base_pieces = gpd.read_file(base_receptors_path)
    for site, zones in base_sites.items():
        temp = base_pieces.loc[np.in1d(base_pieces.zone, zones)]
        temp.loc[:,'grouper'] = 1
        temp = temp.dissolve(by='grouper')
        temp.to_file(os.path.join(outdir,site+'.shp'))

def get_base_sites(): #todo confirm this with MW
    base_sites = {
        # terretorial land authority zone:
        # central is the middle of east and west
        # mid is teh middle of north and south
        'full_tla': range(16),

        'north_tla': [2, 8, 11, 4, 15, 14],
        'north_west_tla': [2, 8],
        'north_central_tla': [4, 11],
        'north_east_tla': [14, 15],

        'mid_tla': [1, 7, 10, 13],
        'mid_west_tla': [1, 7],
        # 'mid_central_tla': [10], # duplicate of urban
        # 'mid_east_tla': [13], # duplicate of urban

        'south_tla': [0, 6, 9, 3, 12, 5],
        'south_west_tla': [0, 6],
        # 'south_central_tla': [9, 3], # duplicate of urban
        # 'south_east_tla': [12, 5], # duplicate of urban

        'west_tla': [2, 8, 1, 7, 6, 0],
        'central_tla': [4, 11, 10, 9, 3],
        'east_tla': [15, 14, 13, 12, 5],

        # urban area:
        'full_urb': [14, 11, 8, 7, 6, 9, 12, 13, 10],

        'west_urb': [8, 7, 6],
        'central_urb': [11, 10, 9],
        'east_urb': [14, 13, 12],

        'north_urb': [8, 11, 14],
        'north_west_urb': [8],
        'north_central_urb': [11],
        'north_east_urb': [14],

        'mid_urb': [7, 10, 13],
        'mid_west_urb': [7],
        'mid_central_urb': [10],
        'mid_east_urb': [13],

        'south_urb': [6, 9, 12],
        'south_west_urb': [6],
        'south_central_urb': [9],
        'south_east_urb': [12]
    }
    return base_sites


def create_interzone_indexes():
    # make a series of loose zones which can be summed togeather to give a full picture of teh zone
    # could save nsmc_num for each one..., which could also minimize the memory in use...
    # return amalgimated sites as well as the indexes to feed into split interzone netcdfs
    base_sites = get_base_sites()
    layer_groups = {'surface': [0],
                    'riccaton': [1],
                    'linwood': [3],
                    'burwood': [5],
                    'wainoni': [7],
                    'shallow': range(4),  # top to bottom of linwood
                    'mid': range(4, 8),  # top of heathcote to bottom of wainoni
                    'deep': range(8, smt.layers),
                    'all_layers': range(smt.layers)
                    }

    sites = {}
    for bs, lg in itertools.product(base_sites.keys(), layer_groups.keys()):
        layers = layer_groups[lg]
        zones = base_sites[bs]
        pieces = ['zone_{}_layer_{}'.format(zone, layer) for zone, layer in itertools.product(zones,layers)]
        sites['{}_{}'.format(lg,bs)] = pieces

    indexes = {}
    base_recpts = smt.shape_file_to_model_array(base_receptors_path, 'zone', True)
    no_flow = smt.get_no_flow().astype(bool)

    for zone in set(base_recpts[np.isfinite(base_recpts)].astype(int)):
        for layer in range(smt.layers):
            t = smt.get_empty_model_grid(True).astype(bool)
            t[layer, np.isclose(base_recpts, zone)] = True
            t[~no_flow] = False
            indexes['zone_{}_layer_{}'.format(zone, layer)] = t

    return indexes, sites


def split_interzone_netcdfs(indir, sites):  # todo
    """

    :param indir: directory with the key netcdfs (e.g. from _get_data from zones)
    :param sites: dictionay {sitename: [component names]}
    :return:
    """
    # aim for somethign similar to split single zone delineation netcdfs, but amalgamate the
    # different views on the zones

    print('splitting netcdf')
    outdir = os.path.join(indir, 'individual_netcdfs')
    if not os.path.exists(os.path.join(indir, 'individual_netcdfs')):
        os.makedirs(outdir)

    data = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        data[name] = nc.Dataset(os.path.join(indir, name + '.nc'))

    for site, components in sites.items():
        outfile = nc.Dataset(os.path.join(outdir, site + '.nc'), 'w')
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
            for suffix in ['', '_cust']:
                temp_var = outfile.createVariable('{}_{}_number{}'.format(name[0:4], name.split('_')[-1][0:2], suffix),
                                                  float,
                                                  ('latitude', 'longitude'),
                                                  fill_value=np.nan)
                temp_var.setncatts({'units': 'number of realisations',
                                    'long_name': '{}_{}{}'.format(name, site, suffix),
                                    'missing_value': np.nan,
                                    'comments': 'number of realisations with particles from given cell'})

                # set up an output array of the correct size (from first component)
                t = data[name].variables['{}{}'.format(components[0], suffix)][:]
                t = t.filled(0).astype(bool)
                t[:] = False
                for comp in components:
                    t = t | (data[name].variables['{}{}'.format(comp, suffix)][:]).filled(0).astype(bool)
                t = t.sum(axis=0).astype(float)
                t[np.isclose(t, 0)] = np.nan
                temp_var[:] = t

        outfile.description = ('source zones for single sources')
        outfile.history = 'created {}'.format(datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.close()


def setup_output_ncs(outdir, sites, model_ids, root_num_part):
    """
    sets up the bulk output netcdfs
    :param outdir: directory to save it all
    :param sites: the site names (e.g. the indexes keys)
    :param model_ids: the model ids this is run for
    :param root_num_part: the cubed
    :return:
    """
    outdata = {}
    for name in ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']:
        outfile = nc.Dataset(os.path.join(outdir, name + '.nc'), 'w')
        x, y = smt.get_model_x_y(False)
        # create dimensions
        outfile.createDimension('latitude', len(y))
        outfile.createDimension('longitude', len(x))
        outfile.createDimension('model_id', len(model_ids))

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

        mid = outfile.createVariable('model_id', str, ('model_id',))
        mid.setncatts({'units': '',
                       'long_name': 'unique_model_identifier',
                       })
        mid[:] = model_ids

        # location initalize the data
        for site in sites:
            temp_var = outfile.createVariable(site, 'u1',
                                              ('model_id', 'latitude', 'longitude'),
                                              fill_value=0)
            temp_var.setncatts({'units': 'bool',
                                'long_name': site,
                                'missing_value': 0,
                                'comments': 'cell in source zone'})

            temp_var = outfile.createVariable('{}_cust'.format(site), 'u1',
                                              ('model_id', 'latitude', 'longitude'),
                                              fill_value=0)
            temp_var.setncatts({'units': 'bool',
                                'long_name': site,
                                'missing_value': 0,
                                'comments': 'cell in source zone'})

        outfile.description = ('source zones for single sources')
        outfile.history = 'created {}'.format(datetime.now().isoformat())
        outfile.source = 'script: {}'.format(sys.argv[0])
        outfile.backward_num_parts = root_num_part ** 3
        outdata[name] = outfile

    return outdata


def add_data_to_nc(out_nc, mid, data, cust_data, run_name):
    mid_idx = np.where(np.array(out_nc.variables['model_id']) == mid)

    # add the non_cust data
    for site in data.keys():
        out_nc.variables[site][mid_idx] = data[site].astype(bool)

    sfr_data, sfr_id_array, unpacked_size, unpacked_shape, losing = cust_data

    # add upstream cust influance
    for site in data.keys():
        temp_sfr_id_array = deepcopy(sfr_id_array)
        temp_sfr_id_array[losing[mid] < 0] = -1

        if not (temp_sfr_id_array[data[site] > 0] >= 0).any():
            out_nc.variables['{}_cust'.format(site)][mid_idx] = data[site].astype(bool)
            continue

        temp_ids = np.array(list(set(temp_sfr_id_array[data[site] > 0]) - {-1})).astype(int)
        # get and unpack the array
        temp = np.unpackbits(sfr_data[run_name][mid])[:unpacked_size].reshape(unpacked_shape).astype(bool)
        temp = temp[:temp_ids.max() + 1].sum(axis=0)
        temp = temp.astype(bool) | data[site].astype(bool)
        out_nc.variables['{}_cust'.format(site)][mid_idx] = temp


def _get_data_for_zones(outdata, run_name, model_ids, indexes, root_num_part, recalc_backward_tracking):
    """
    get and amalgamate up the data for the zone delination
    :param outdata: dictionary of keys ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']
                    and values open netcdf file to write data in
    :param run_name: the name to call this run of backward models (to prevent overwrite)
    :param model_ids: list of model ids
    :param indexes: dictionary of boolean arrays for the targets
    :param root_num_part: the cubic root of the number of particles to release in each backward modpath cell
    :param recalc_backward_tracking: bool, if true re-run the backward partical tracking
    :return: amalg_weak_forward, amalg_strong_forward, amalg_strong_back, amalg_weak_back, forward_strongs_num_parts
    """
    assert {'forward_weak', 'forward_strong', 'backward_strong', 'backward_weak'} == set(outdata.keys())
    assert all([isinstance(e, nc.Dataset) for e in outdata.values()])
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = os.path.join(get_base_results_dir('backward', socket.gethostname()), run_name)

    cust_data = get_cust_mapping(run_name, model_ids)

    # backward weak
    print('calculating backward weak\n\n')
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i + 1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'weak', mid),
                                           mp_name='{}_weak'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=True,
                                           recalc=recalc_backward_tracking)
        add_data_to_nc(outdata['backward_weak'], mid, temp, cust_data, 'backward_weak')

    # backward strong
    print('calculating backward strong\n\n')
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i + 1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'strong', mid),
                                           mp_name='{}_strong'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=False,
                                           recalc=recalc_backward_tracking)
        add_data_to_nc(outdata['backward_strong'], mid, temp, cust_data, 'backward_strong')

    # forward weak
    print('calculating forward weak\n\n')
    f_em_paths = get_forward_emulator_paths(model_ids, True)

    for i, (mid, path) in enumerate(f_em_paths.items()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i + 1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        add_data_to_nc(outdata['forward_weak'], mid, temp, cust_data, 'forward_weak')

    # forward strong
    print('calculating forward strong\n\n')
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)

    for i, (mid, path) in enumerate(f_em_paths.items()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i + 1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        add_data_to_nc(outdata['forward_strong'], mid, temp, cust_data, 'forward_strong')


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

    # save the data
    outdata = setup_output_ncs(outdir=outdir,
                               model_ids=model_ids,
                               root_num_part=root_num_part, sites=indexes.keys())

    _get_data_for_zones(outdata, run_name, model_ids, indexes,
                        root_num_part, recalc_backward_tracking)

    return outdata


def run_interzone_source_zones(recalc=False, recalc_backward_tracking=False):
    """
    run the singel source zones... final wrapper
    :param recalc: bool if True recalculate the netcdf file
    :param recalc_backward_tracking: bool if True recalculate the modpath particle tracking simulations
    :return:
    """
    indexes, sites = create_interzone_indexes()
    base_outdir = r"C:\mh_waimak_models\interzone_source_zones"
    print('running for AshOpt')
    outdir = os.path.join(base_outdir, 'AshOpt')
    create_zones(model_ids=['AshOpt'], run_name='AshOpt_interzone_sources',
                 outdir=outdir, root_num_part=1,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_interzone_netcdfs(outdir, sites)

    print('running for 165 models')
    outdir = os.path.join(base_outdir, 'stocastic set')
    stocastic_model_ids = get_stocastic_set()
    create_zones(model_ids=stocastic_model_ids,
                 run_name='stocastic_set_interzone_sources',
                 outdir=outdir, root_num_part=1,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_interzone_netcdfs(outdir, sites)


if __name__ == '__main__':
    run_interzone_source_zones(recalc=True,recalc_backward_tracking=False)
    print('done')
