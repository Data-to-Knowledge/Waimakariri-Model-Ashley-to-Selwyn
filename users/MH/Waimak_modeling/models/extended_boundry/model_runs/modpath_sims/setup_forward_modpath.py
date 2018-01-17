# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/12/2017 1:47 PM
"""

from __future__ import division
from core import env
import numpy as np
import flopy
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.modpath_wrapper import \
    create_mp_slf, get_cbc
import os
import pandas as pd


def part_group_cell_mapper(bd_type):
    js, iss = np.meshgrid(range(smt.cols), range(smt.rows))  # zero indexed to agree with python interpretation
    idx = bd_type.flatten() != -1
    out = dict(zip(range(1, idx.sum() + 1), list(zip(iss.flatten()[idx], js.flatten()[idx]))))
    return out


def make_mp_forward_particles(cbc_path, min_part=1, max_part=None):
    """
    make modpath particle locations from the cbc file.  particles are created in each top layer cell with the number
    relative to the influx
    :param cbc_path: path to the cell by cell budget file
    :param min_part: the minimum number of particles in each cell
    :param max_part: the number above which particles will be truncated or None (no truncation)
    :return: np.record array
    """
    # particles must have a label as a string otherwise it terminates
    # mass weighted particles for each of the models for the sfr/well(races/rivers) and rch

    print('reading cbc')
    # calculate the flow into each top layer cell
    indata = flopy.utils.CellBudgetFile(cbc_path)
    rch = indata.get_data(kstpkper=(0, 0), text='recharge', full3D=True)[0][0].filled(0)  # take only the first row
    # take only the first row, no particles in southwertern boundary for wells
    well = indata.get_data(kstpkper=(0, 0), text='wells', full3D=True)[0][0].filled(0)
    sfr = indata.get_data(kstpkper=(0, 0), text='stream leakage', full3D=True)[0][0].filled(
        0)  # take only the first row
    well[well < 0] = 0
    sfr[sfr < 0] = 0
    flow = rch + well + sfr
    flow[smt.get_no_flow(0) != 1] = 0

    # generate particles (minimum of 1 per cell)
    num_parts = np.round(flow / flow[flow > 0].min() * min_part).astype(int)
    if max_part is not None:
        num_parts[num_parts > max_part] = int(max_part)

    # identify boundary condition types
    bd_type = smt.get_empty_model_grid()  # 0=rch,1=well,2==sfr
    temp = np.concatenate((rch[np.newaxis], well[np.newaxis], sfr[np.newaxis]), axis=0).max(axis=0)
    bd_type[np.isclose(rch, temp)] = 0
    bd_type[np.isclose(well, temp)] = 1
    bd_type[np.isclose(sfr, temp)] = 2
    bd_type[(smt.get_no_flow(0) != 1) | np.isclose(temp, 0)] = -1

    outdata = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(num_parts.sum())
    outdata['label'] = 's'
    outdata['k0'] = 0  # I think that particles in flopy are 1 indexed, but 0 means active most layer
    js, iss = np.meshgrid(range(1, smt.cols + 1), range(1, smt.rows + 1))  # this is passed to the file as 1 indexed

    idx = bd_type.flatten() != -1
    group_dict = part_group_cell_mapper(bd_type)
    start_idx = 0
    print('generating particles')
    for l, (num, i, j, bt) in enumerate(
            zip(num_parts.flatten()[idx], iss.flatten()[idx], js.flatten()[idx], bd_type.flatten()[idx])):
        if num == 0:
            raise ValueError('unexpected zero points')
        end_idx = start_idx + num
        # set particle starting location
        outdata['particlegroup'][start_idx:end_idx] = l + 1
        outdata['groupname'][start_idx:end_idx] = '{:03d}_{:03d}'.format(*group_dict[l + 1])
        outdata['i0'][start_idx:end_idx] = i
        outdata['j0'][start_idx:end_idx] = j
        outdata['xloc0'][start_idx:end_idx] = np.random.uniform(size=num)
        outdata['yloc0'][start_idx:end_idx] = np.random.uniform(size=num)
        outdata['zloc0'][start_idx:end_idx] = 1 if bt == 0 else np.random.uniform(size=num)
        start_idx = end_idx
    return outdata, bd_type, num_parts


def setup_run_forward_modpath(cbc_path, mp_ws, mp_name, min_part=1, max_part=None, capt_weak_s=False):
    """
    a wrapper to quickly setup and run a forward modpath simulation for all influx cells
    :param cbc_path: the path to the cell budget file, it is expected that the same pattern is used for the heads file
    :param mp_ws: the directory to put the model path into
    :param mp_name: the name of the model the model id is not passed to the name here
    :param min_part: the minimum number of particles in a cell
    :param max_part: the number above which particles will be truncated or None (no truncation)
    :param capt_weak_s: Bool if True terminate particles in weak source/sinks
    :return:
    """
    if not os.path.exists(mp_ws):
        os.makedirs(mp_ws)
    particles, bd_type, num_parts = make_mp_forward_particles(cbc_path, min_part=min_part, max_part=max_part)
    particles = pd.DataFrame(particles)
    np.savetxt(os.path.join(mp_ws, '{}_bnd_type.txt'.format(mp_name)), bd_type)
    np.savetxt(os.path.join(mp_ws, '{}_num_parts.txt'.format(mp_name)), num_parts)

    temp_particles = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(0)
    mp = create_mp_slf(particle_data=temp_particles, mp_ws=mp_ws, hdfile=cbc_path.replace('cbc', 'hds'),
                       budfile=cbc_path, disfile=cbc_path.replace('cbc', 'dis'), mp_name=mp_name,
                       capt_weak_s=capt_weak_s)
    print('writing model {}; ignore the "no data to write" comment (this is a hack)'.format(mp_name))
    mp.write_input()
    mp.write_name_file()

    # write loc file with pandas to save time
    # simple speed test writing particles with flopy and running model took 30 min, writing with pandas took __min
    loc_path = mp.get_package('loc').fn_path
    loc_package = mp.get_package('loc')
    # write groups
    print ('writing group loc data')
    groups = particles[['particlegroup', 'groupname']].groupby('particlegroup').count().reset_index().rename(
        columns={'groupname': 'count'})
    groups.loc[:, 'groupname'] = groups.loc[:, 'particlegroup'].replace(
        dict(particles[['particlegroup', 'groupname']].itertuples(False, None)))
    group_count = len(groups.index)
    groups = pd.Series(groups[['groupname', 'count']].astype(str).values.flatten())
    with open(loc_path, 'w') as f:
        f.write('{}\n'.format(loc_package.heading))
        f.write('{:d}\n'.format(loc_package.input_style))
        f.write('{}\n'.format(group_count))

    groups.to_csv(loc_path, False, sep=' ', header=False, mode='a')

    # write particle data
    print ('writing loc particle data')
    particles.drop('groupname', 1, inplace=True)
    particles.to_csv(loc_path, sep=' ', header=False, index=False, mode='a')

    print('running model {}'.format(mp_name))
    mp.run_model()


if __name__ == '__main__':
    import time

    t = time.time()
    test = setup_run_forward_modpath(r"C:\Users\MattH\Desktop\NsmcBase_simple_modpath\NsmcBase_modpath_base.cbc",
                                     r"C:\Users\MattH\Desktop\test_write\part2", 'one_part_per_cell')
    print(time.time() - t)
    print('done')
