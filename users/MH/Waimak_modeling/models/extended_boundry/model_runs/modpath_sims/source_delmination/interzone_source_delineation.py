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
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import get_stocastic_set
from single_zone_delination import get_modeflow_dir_for_source, get_base_results_dir, get_cust_mapping, define_source_from_backward, get_cbc, get_forward_emulator_paths, define_source_from_forward


def create_interzone_indexes():
    # make a series of loose zones which can be summed togeather to give a full picture of teh zone
    #todo shit these can't be summed togeather due to overlaps
    # could save nsmc_num for each one..., which could also minimize the memory in use...

    raise NotImplementedError

def split_interzone_netcdfs(outdir):
    # aim for somethign similar to split single zone delineation netcdfs, but amalgamate the
    #different views on the zones
    raise NotImplementedError

def setup_output_ncs(outdir,model_ids,root_num_part):
    versions = ['forward_weak', 'forward_strong', 'backward_strong', 'backward_weak']

    #todo setup the 4 netcdfs, which will hold the model_id number (or similar or boolean)
    raise NotImplementedError
    return outdata

def add_data_to_nc(out_nc, mid, data, cust_data):
    #todo add the data to the netcdf includign the cust data
    # todo how to implement cust!!!
    #todo I may have to rehash so that out_nc goes to outdata (e.g. dict of ncs)
    raise NotImplementedError

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
    assert all([isinstance(e,nc.Dataset) for e in outdata.values()])
    modflow_dir = get_modeflow_dir_for_source()
    backward_dir = os.path.join(get_base_results_dir('backward', socket.gethostname()), run_name)

    cust_data = get_cust_mapping(run_name, model_ids)

    # backward weak
    print('calculating backward weak\n\n')
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i+1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'weak', mid),
                                           mp_name='{}_weak'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=True,
                                           recalc=recalc_backward_tracking)
        add_data_to_nc(outdata['backward_weak'], mid, temp, cust_data)

    # backward strong
    print('calculating backward strong\n\n')
    for i, mid in enumerate(model_ids):
        print('model: {}, {} of {}'.format(mid, i+1, len(model_ids)))
        temp = define_source_from_backward(indexes,
                                           mp_ws=os.path.join(backward_dir, 'strong', mid),
                                           mp_name='{}_strong'.format(mid),
                                           cbc_file=get_cbc(model_id=mid, base_dir=modflow_dir),
                                           root3_num_part=root_num_part, capt_weak_s=False,
                                           recalc=recalc_backward_tracking)
        add_data_to_nc(outdata['backward_strong'], mid, temp, cust_data)

    # forward weak
    print('calculating forward weak\n\n')
    f_em_paths = get_forward_emulator_paths(model_ids, True)

    for i, (mid, path) in enumerate(f_em_paths.items()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i+1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        add_data_to_nc(outdata['forward_weak'], mid, temp, cust_data)

    # forward strong
    print('calculating forward strong\n\n')
    f_em_paths = get_forward_emulator_paths(model_ids, weak_sink=False)

    for i, (mid, path) in enumerate(f_em_paths.items()):
        print('{}, {} of {}'.format(os.path.basename(path[0]), i+1, len(f_em_paths)))
        temp = define_source_from_forward(emulator_path=path[0], bd_type_path=path[1], indexes=indexes)
        add_data_to_nc(outdata['forward_weak'], mid, temp, cust_data)


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
                             root_num_part=root_num_part)

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
    indexes = create_interzone_indexes()
    base_outdir = r"C:\mh_waimak_models\interzone_source_zones"
    print('running for AshOpt')
    outdir = os.path.join(base_outdir, 'AshOpt')
    create_zones(model_ids=['AshOpt'], run_name='AshOpt_interzone_sources',
                 outdir=outdir, root_num_part=1,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_interzone_netcdfs(outdir)

    print('running for 165 models')
    outdir = os.path.join(base_outdir, 'stocastic set')
    stocastic_model_ids = get_stocastic_set()
    create_zones(model_ids=stocastic_model_ids,
                 run_name='stocastic_set_interzone_sources',
                 outdir=outdir, root_num_part=1,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_interzone_netcdfs(outdir)
