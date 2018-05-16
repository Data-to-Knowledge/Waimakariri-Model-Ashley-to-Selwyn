# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/04/2018 1:41 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.modpath_wrapper import \
    create_mp_slf, get_cbc, export_paths_to_shapefile
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.setup_reverse_modpath import particle_loc_from_grid
import pandas as pd
from glob import glob
import os
import numpy as np
import flopy

if __name__ == '__main__':
    temp = smt.get_empty_model_grid(True)
    temp[0] = np.isfinite(smt.shape_file_to_model_array(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\m_ex_bd_inputs\raw_sw_samp_points\sfr\swaz\eyre_swaz.shp", 'stop',True))
    grid_locs = smt.model_where(temp)
    groups = np.repeat(['1'],len(grid_locs))
    group_mapper, particles = particle_loc_from_grid(grid_locs, groups, root3_num_part=10)
    outdir ="C:\Users\MattH\Desktop\eyre_modpath"
    cbc = get_cbc('NsmcBase',outdir)
    temp_particles = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(0)
    mp = create_mp_slf(temp_particles, mp_ws=outdir, hdfile=cbc.replace('.cbc','.hds'), budfile=cbc,
                       disfile=cbc.replace('.cbc','.dis'),
                    prsity=0.3, prsitycb=0.3, mp_name='eyre_river')
    mp.write_input()
    print('writing model {}; ignore the "no data to write" comment (this is a hack)')
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


    mp.run_model()
    path_file = glob(os.path.join(outdir, '*.mppth'))[0]
    export_paths_to_shapefile(path_file,
                              os.path.join(outdir,'check_eyre.shp'))
