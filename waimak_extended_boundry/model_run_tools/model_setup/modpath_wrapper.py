# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/12/2017 3:06 PM
"""

from __future__ import division
import flopy_mh as flopy
import numpy as np
from env import sdp
import os
from waimak_extended_boundry import smt
from base_modflow_wrapper import import_gns_model
from waimak_extended_boundry.model_run_tools.convergance_check import \
    modflow_converged
from traceback import format_exc


# make something to run modpath simulations with particles in the top most active cells.


def get_cbc_mp(kwargs):
    try:
        get_cbc(**kwargs)
        success = 'converged'
    except Exception:
        success = format_exc()
    return kwargs['model_id'], success


def get_cbc(model_id, base_dir, recalc=False):
    # if I wanted to split layers I would need to do it here...
    cbc_path = os.path.join(base_dir, '{}_for_modpath'.format(model_id), '{}_for_modpath.cbc'.format(model_id))

    if os.path.exists(cbc_path) and not recalc:
        return cbc_path

    m = import_gns_model(model_id, 'for_modpath', os.path.join(base_dir, 'for_modpath'), False, set_hdry=True)
    m.write_name_file()
    m.write_input()
    m.run_model()
    con = modflow_converged(os.path.join(m.model_ws, '{}.list'.format(m.name)))
    if not con:
        raise ValueError('model did not converge')
    return cbc_path


def create_mp_slf(particle_data, m=None, mp_ws=None, hdfile=None, budfile=None, disfile=None,
                  prsity=0.3, prsitycb=0.3, mp_name=None, direction='forward',
                  simulation_type='pathline', capt_weak_s=False, time_pts=1, hnoflo=1e+30, hdry=-888.0,
                  laytype=(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)):
    """
    create a modpath simulation which derives particles from point data.  assumptions are listed below
    :param particle_data: record array of particle data, Note that the label must be explicitly passed (cannot be '')
    :param m: modflow object which has been run or None
    :param mp_ws: the directory to save the modpath files and output
    :param hdfile: the full path to the head file, or None
    :param budfile: the full path to the budget file or None
    :param disfile: the full path to the dis file or None
    :param prsity: porosity
    :param prsitycb: porosity of the confining layer
    :param mp_name: the name for the modpath model if None then mirrors the modflow model name
    :param direction: tracking direction 'forward' or 'backward'
    :param simulation_type: one of [endpoint, pathline, timeseries]
    :param capt_weak_s: if True weak sources/sinks capture particles
    :param time_pts:  either int of number of time points to calculate or array of time points see TimePointOption for
                     documentation
    :return: modpath model
    """
    model_passed = (isinstance(m, flopy.modflow.Modflow) and hdfile is None and budfile is None and disfile is None)
    strs_passed = (m is None and all([isinstance(e, str) for e in [mp_ws, mp_name, hdfile, budfile, disfile, ]]))
    assert model_passed or strs_passed, 'either the model or paths must be defined, see documentation'

    if m is not None:
        hdiu = m.oc.iuhead
        hdfile = os.path.join(m.model_ws, m.get_output(unit=hdiu))
        budfile = os.path.join(m.model_ws, m.get_output(unit=hdiu).replace('.hds', '.cbc'))
        disfile = os.path.join(m.model_ws, m.get_package('dis').file_name[0])

        if mp_ws is None:
            mp_ws = m.model_ws

        if mp_name is None:
            mp_name = '{}_mp'.format(m.name)
    else:
        m = flopy.modflow.Modflow(version='mfnwt')
        flopy.modflow.ModflowDis(m, smt.layers, smt.rows, smt.cols, 1)
        flopy.modflow.ModflowUpw(m)

    mp = flopy.modpath.Modpath(modelname=mp_name,
                               exe_name="{}/models_exes/modpath.6_0/bin/mp6.exe".format(sdp),
                               modflowmodel=m,
                               model_ws=mp_ws,
                               listunit=6,
                               dis_file=disfile,
                               head_file=hdfile,
                               budget_file=budfile)

    mpb = flopy.modpath.ModpathBas(mp,
                                   hnoflo=hnoflo,
                                   hdry=hdry,
                                   bud_label=None,  # what should this be? might be the header it seems to work as None
                                   laytyp=laytype,
                                   ibound=smt.get_no_flow(),
                                   prsity=prsity,
                                   prsityCB=prsitycb
                                   )

    # pass argments to sim codes
    if simulation_type == 'endpoint':
        sim_type = 1
    elif simulation_type == 'pathline':
        sim_type = 2
    elif simulation_type == 'timeseries':
        sim_type = 3
    else:
        raise ValueError('unexpected simulation type: {}'.format(simulation_type))

    if direction == 'forward':
        t_dir = 1
    elif direction == 'backward':
        t_dir = 2
    else:
        raise ValueError('unexpected tracking direction: {}'.format(direction))

    if capt_weak_s:
        wk_s = 2
    else:
        wk_s = 1

    if simulation_type == 'endpoint':
        time_pt_o = 1
        time_ct = 1
    else:
        if isinstance(time_pts, int):
            time_pt_o = 2
            time_ct = time_pts
        elif time_pts is None:
            raise ValueError('timepoints cannot be None')
        else:
            time_pts = np.atleast_1d(time_pts)
            time_pt_o = 3
            time_ct = len(time_pts)
            time_pts = list(time_pts)

    sim = flopy.modpath.ModpathSim(mp,
                                   option_flags=[sim_type,  # simulation type
                                                 t_dir,  # tracking direction
                                                 wk_s,  # weak sink option
                                                 wk_s,  # weak source option
                                                 2,  # reference time option (this is how the previous term speficed)
                                                 2,  # StopOption extend model to termination points
                                                 2,  # read particles from file
                                                 time_pt_o,  # time point options as above
                                                 1,  # default budget option from create sim
                                                 1,  # set all zones to 1
                                                 1,  # no retardation factors are read
                                                 1],  # advective observations not computed or saved
                                   ref_time=0,  # default
                                   ref_time_per_stp=[0, 0, 1.0],  # default
                                   stop_time=None,  # default
                                   group_name=['group_1'],  # default
                                   group_placement=[[1,  # must be 1
                                                     1,  # not used if reading particles from file
                                                     1,  # placement option
                                                     0,  # relaese time
                                                     1,  # release option
                                                     1]],  # particles not generated for constant head boundarys
                                   release_times=[[1, 1]],  # not using as releasing all at once
                                   group_region=[[1, 1, 1, 1, 1, 1]],  # not using as reading particles from file
                                   mask_nlay=[1],  # not using as reading particles from file
                                   mask_layer=[1],  # not using as reading particles from file
                                   mask_1lay=[1],  # probably not using
                                   face_ct=[1],  # not using as reading particles from file
                                   ifaces=[[6, 1, 1]],  # default
                                   part_ct=[[1, 1, 1]],  # not using as reading particles from file
                                   time_ct=time_ct,
                                   release_time_incr=1,  # default should calculate in daily increments
                                   time_pts=time_pts,
                                   particle_cell_cnt=[[2, 2, 2]],  # not using as reading particles from file
                                   cell_bd_ct=1,  # default not used
                                   bud_loc=[[1, 1, 1, 1]],  # default not used
                                   trace_id=1,  # default not used
                                   stop_zone=1,  # default not used
                                   zone=1,  # default
                                   retard_fac=1.0,  # default
                                   retard_fcCB=1.0)  # default

    start_loc = flopy.modpath.mpsim.StartingLocationsFile(mp)
    start_loc.data = particle_data

    return mp


def export_paths_to_shapefile(paths_file, shape_file, particle_ids=None):
    """
    export paths to a shapefile
    :param paths_file:
    :param shape_file:
    :param particle_ids:
    :return:
    """
    # generate spatial reference
    spatial_ref = flopy.utils.SpatialReference(delr=np.full((smt.rows), 200), delc=np.full((smt.cols), 200), lenuni=2,
                                               xul=smt.ulx, yul=smt.uly, rotation=smt.rotation,
                                               proj4_str="EPSG:2193", units='meters',
                                               length_multiplier=1.)
    paths = flopy.utils.PathlineFile(paths_file)
    if particle_ids is None:
        pathdata = [e[['x', 'y', 'z', 'k', 'id']] for e in paths.get_alldata()]
        paths.write_shapefile(shpname=shape_file, sr=spatial_ref)
    else:
        particle_ids = np.atleast_1d(particle_ids)
        pathdata = [paths.get_data(e) for e in particle_ids]
        paths.write_shapefile(particle_data=pathdata, shpname=shape_file, sr=spatial_ref)

def pathline_file_to_hdf():
    # tod this would cut the stored size in half with no compression...  plus you can query data within the store
    raise NotImplementedError


if __name__ == '__main__':
    run_model = False
    if run_model:
        particle_data = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(2)
        particle_data['particleid'] = [120, 150]
        particle_data['i0'][:] = [120, 150]
        particle_data['j0'][:] = [120, 150]
        particle_data['particlegroup'][:] = 1
        particle_data['groupname'][:] = 'part1'
        particle_data['label'][:] = 'test'
        model_path = r"C:\Users\MattH\Desktop\NsmcBase_modpath_tester\NsmcBase_modpath_tester.nam"
        m = flopy.modflow.Modflow.load(model_path, model_ws=os.path.dirname(model_path))
        mp = create_mp_slf(particle_data, m)
        mp.write_input()
        mp.write_name_file()
        mp.run_model()
    export_paths_to_shapefile(r"C:\Users\MattH\Desktop\NsmcBase_modpath_tester\NsmcBase_modpath_tester_mp.mppth",
                              r"C:\Users\MattH\Desktop\NsmcBase_modpath_tester\test_shape.shp", particle_ids=[119])
    print('done')
