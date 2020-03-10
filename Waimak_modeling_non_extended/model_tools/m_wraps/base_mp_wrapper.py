"""
Author: matth
Date Created: 20/04/2017 11:08 AM
"""

from __future__ import division
import numpy as np
import itertools
import flopy_mh as flopy
from warnings import warn
from users.MH.Waimak_modeling.supporting_data_path import sdp
from users.MH.Waimak_modeling.model_tools.basic_tools import calc_elv_db, convert_coords_to_matix


def create_mp(m, start_time, trackdir='forward', packages='RCH', simtype='pathline', time_pts = None, capt_weak_s=False,
              prsity=0.3, prsityCB = 0.3, ParticleColumnCount=4, ParticleRowCount=4, mp_name = None):
    """
    wrapper to set defaults for many of the modflow options I need pass in a modflow model which has run
    :param m: modflow object which has been run
    :param start_time: Sets the value of MODPATH reference time relative to MODFLOW time.
                       float: value of MODFLOW simulation time at which to start the particle tracking simulation.
                              Sets the value of MODPATH ReferenceTimeOption to 1.
                       tuple: (period, step, time fraction) MODFLOW stress period, time step and fraction
                              between 0 and 1 at which to start the particle tracking simulation.
                              Sets the value of MODPATH ReferenceTimeOption to 2.
    :param trackdir: tracking direction forward or reverse
    :param packages: package to initiate particles at supported: RCH, WEL, MNW2
    :param simtype: Keyword defining the MODPATH simulation type. Available: 'endpoint' 'pathline' and 'timeseries'
    :param time_pts: either int of number of time points to calculate or array of time points see TimePointOption for
                     documentation
    :param capt_weak_s: if True weak sources/sinks capture particles
    :param prsity: porosity
    :param prsityCB: porosity of the confining layer
    :param ParticleColumnCount: Columns of particles to start on each cell index face (iface)
    :param ParticleRowCount: Rows of particles to start on each cell index face (iface).
    :param mp_name: the name for the modpath model if None then mirrors the modflow model name
    :return: modpath model with established simulation
    """

    if mp_name is None:
        mp_name = '{}_mp'.format(m.name)

    mp = flopy.modpath.Modpath(modelname=mp_name,
                               exe_name="{}/models_exes/modpath.6_0/bin/mp6.exe".format(sdp),
                               modflowmodel=m,
                               model_ws=m.model_ws,
                               listunit=6
                               )

    mpb = flopy.modpath.ModpathBas(mp,
                                   hnoflo=m.bas6.hnoflo,
                                   hdry=m.lpf.hdry,
                                   bud_label=None,  # what should this be? might be the header it seems to work as None
                                   laytyp=m.lpf.laytyp.array,
                                   ibound=m.bas6.ibound.array,
                                   prsity=prsity,
                                   prsityCB=prsityCB
                                   )

    # 1 indexed options (e.g. the first entry in the options list is 1 here)
    # 'SimulationType': 1,
    # 'TrackingDirection': 2,
    # 'WeakSinkOption': 3,
    # 'WeakSourceOption': 4,
    # 'ReferenceTimeOption': 5,
    # 'StopOption': 6,
    # 'ParticleGenerationOption': 7, # this is how you read from paricle start file
    # 'TimePointOption': 8, #less sure about this one
    # 'BudgetOutputOption': 9,
    # 'ZoneArrayOption': 10}
    # 'RetardationOption': 11,
    # 'AdvectiveObservationsOption': 12,

    sim = mp.create_mpsim(trackdir=trackdir,
                          simtype=simtype,
                          packages=packages,
                          start_time=start_time,
                          ParticleColumnCount=ParticleColumnCount,
                          ParticleRowCount=ParticleRowCount,
                          )
    # set weak sink/source options
    if capt_weak_s:
        sim.options_dict['WeakSinkOption'] = 2
        sim.options_dict['WeakSourceOption'] = 2
        sim.option_flags[2:4] = [2,2]
    else:
        sim.options_dict['WeakSinkOption'] = 1
        sim.options_dict['WeakSourceOption'] = 1
        sim.option_flags[2:4] = [1,1]

    # set timepoint values # check this
    if simtype == 'endpoint':
        sim.options_dict['TimePointOption'] = 1
        sim.option_flags[7] = 1
    else:
        if isinstance(time_pts, int):
            sim.options_dict['TimePointOption'] = 2
            sim.option_flags[7] = 2
            sim.time_ct = time_pts
        elif time_pts is None:
            warn('no time points passed, will use defaults from flopy.modpath.Modpath.create_mpsim()')
        else:
            time_pts = np.atleast_1d(time_pts)
            sim.options_dict['TimePointOption'] = 3
            sim.option_flags[7] = 3
            sim.time_ct = len(time_pts)
            sim.time_pts = list(time_pts)

    return mp

def create_mp_slf (m, particle_data, prsity=0.3, prsityCB = 0.3, mp_name=None, direction='forward',
                   simulation_type='endpoint', capt_weak_s=False, time_pts=1): # debug this
    """
    create a modpath simulation which derives particles from point data.  assumptions are listed below
    :param m: modflow object which has been run
    :param particle_data: record array of particle data
    :param prsity: porosity
    :param prsityCB: porosity of the confining layer
    :param mp_name: the name for the modpath model if None then mirrors the modflow model name
    :param direction: tracking direction forward or backward
    :param simulation_type: one of [endpoint, pathline, timeseries]
    :param capt_weak_s: if True weak sources/sinks capture particles
    :param time_pts:  either int of number of time points to calculate or array of time points see TimePointOption for
                     documentation
    :return: modpath model
    """
    if mp_name is None:
        mp_name = '{}_mp'.format(m.name)

    mp = flopy.modpath.Modpath(modelname=mp_name,
                               exe_name="{}/models_exes/modpath.6_0/bin/mp6.exe".format(sdp),
                               modflowmodel=m,
                               model_ws=m.model_ws,
                               listunit=6
                               )

    mpb = flopy.modpath.ModpathBas(mp,
                                   hnoflo=m.bas6.hnoflo,
                                   hdry=m.lpf.hdry,
                                   bud_label=None,  # what should this be? might be the header it seems to work as None
                                   laytyp=m.lpf.laytyp.array,
                                   ibound=m.bas6.ibound.array,
                                   prsity=prsity,
                                   prsityCB=prsityCB
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
        t_dir =2
    else:
        raise ValueError('unexpected tracking direction: {}'.format(direction))

    if capt_weak_s:
        wk_s = 2
    else:
        wk_s =1

    if simulation_type == 'endpoint':
        time_pt_o = 1
        time_ct = 1 # this may cause problems
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
                                   option_flags=[sim_type, # check these are right against dictionary
                                                 t_dir,
                                                 wk_s,
                                                 wk_s,
                                                 2,  # this is how the previous term speficed referenc option
                                                 2,  # extend model as far forward through termination points
                                                 2,  # read particles from file
                                                 time_pt_o, # time point options as above
                                                 1,  # default budget option from create sim
                                                 1,  # set all zones to 1
                                                 1,  # no retardation factors are read
                                                 1],# advective obeservations not computed or saved
                                   ref_time=0,# default
                                   ref_time_per_stp=[0, 0, 1.0],# default
                                   stop_time=None,# default
                                   group_name=['group_1'],# default
                                   group_placement=[[1, # must be 1
                                                     1, # not used if reading particles from file
                                                     1, # placement option
                                                     0, # relaese time
                                                     1, #  release option
                                                     1]], # particles not generated for constant head boundarys
                                   release_times=[[1, 1]], # not using as releasing all at once
                                   group_region=[[1, 1, 1, 1, 1, 1]],# not using as reading particles from file
                                   mask_nlay=[1],# not using as reading particles from file
                                   mask_layer=[1],# not using as reading particles from file
                                   mask_1lay=[1],# probably not using
                                   face_ct=[1], # not using as reading particles from file
                                   ifaces=[[6, 1, 1]],#default
                                   part_ct=[[1, 1, 1]], # not using as reading particles from file
                                   time_ct=time_ct,
                                   release_time_incr=1, #default should calculate in daily increments
                                   time_pts=time_pts,
                                   particle_cell_cnt=[[2, 2, 2]],# not using as reading particles from file
                                   cell_bd_ct=1,#default not used
                                   bud_loc=[[1, 1, 1, 1]],#default not used
                                   trace_id=1,# default not used
                                   stop_zone=1,# default not used
                                   zone=1,#default
                                   retard_fac=1.0,# default
                                   retard_fcCB=1.0)#default

    start_loc = flopy.modpath.mpsim.StartingLocationsFile(mp)
    start_loc.data = particle_data

    return mp

def particle_loc_from_grid(grid_locs, root3_num_part = 2):
    """

    :param grid_locs: a list of tuples(k,i,j) or tuple (k,i,j)
    :return: starting location data
    """
    # grid_locs a list of tuples(k,i,j)
    # take a list of cell locations (k,i,j) and return a record array for starting location data
    # assume 8 particles per cel for now
    grid_locs = list(np.atleast_1d(grid_locs))
    num_per_cel = root3_num_part**3
    num_par = len(grid_locs)*num_per_cel
    locals_vals = [[e/(root3_num_part+1)] for e in range(1,root3_num_part+1)]
    ploc = list(itertools.product(grid_locs,locals_vals,locals_vals,locals_vals))
    ploc = np.array([list(itertools.chain(*e)) for e in ploc])
    # k,i,j,x,y,z
    outdata = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(npt=num_par)

    outdata['k0'] = ploc[:, 0]
    outdata['i0'] = ploc[:, 1]
    outdata['j0'] = ploc[:, 2]
    outdata['xloc0'] = ploc[:, 3]
    outdata['yloc0'] = ploc[:, 4]
    outdata['zloc0'] = ploc[:, 5]
    outdata['initialtime'] = 1 # should not matter for now as I am releasing all at once
    # for future note that this needs a lable string
    return outdata

def particle_ring_around_screens(radius, n_pts, xs, ys, screen_zs, labels=None, forgive=False):
    circle_points = np.array([
        (radius * np.cos(theta), radius * np.sin(theta))
        for theta in (np.pi*2 * i/n_pts for i in range(n_pts))
    ])

    if labels is None:
        labels = np.zeros(np.array(xs).shape,dtype='S1')
        labels[:] = ''
    nparticles = len(list(sum(screen_zs, ())))*n_pts

    outdata = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(npt=nparticles)
    elv_db = calc_elv_db()
    i=0
    for x, y, zs, lab in zip(xs, ys, screen_zs, labels):
        temp_p_points = circle_points + np.array([x, y])[np.newaxis,:]

        for z in zs:
            for pt in temp_p_points:
                if forgive:
                    try:
                        zloc, yloc, xloc = convert_coords_to_matix(pt[0], pt[1], z, elv_db=elv_db, return_AE=True)
                    except ValueError as val:
                        continue
                else:
                    zloc, yloc, xloc = convert_coords_to_matix(pt[0], pt[1], z, elv_db=elv_db, return_AE=True)


                outdata['k0'][i] = zloc[0]
                outdata['i0'][i] = yloc[0]
                outdata['j0'][i] = xloc[0]
                outdata['zloc0'][i] = zloc[1]
                outdata['yloc0'][i] = yloc[1]
                outdata['xloc0'][i] = xloc[1]
                outdata['label'][i] = lab
                outdata['groupname'][i] = 'g1'
                i += 1

    outdata = outdata[outdata['groupname']=='g1']


    outdata['initialtime'] = 1  # should not matter for now as I am releasing all at once

    return outdata

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from users.MH.Waimak_modeling.model_tools.well_values import get_all_well_data
    well_data = get_all_well_data().head(1)

    temp = particle_ring_around_screens(25, 5, well_data['x'],well_data['y'], well_data['screen_z'])

    plt.scatter(temp['j0'] + temp['xloc0'], -1*(temp['i0']-temp['yloc0']))
    plt.axes().set_aspect('equal', 'datalim')

    grid_pts = [(well_data.loc[e,'layer'],well_data.loc[e,'row'],well_data.loc[e,'col']) for e in well_data.index]

    temp2 = particle_loc_from_grid(grid_pts)
    print('done')
