# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/04/2018 1:07 PM
"""

from __future__ import division
import os
import numpy as np
import flopy_mh as flopy

from waimak_extended_boundry import smt
from waimak_extended_boundry import mod_gns_model
from waimak_extended_boundry import get_base_well
from waimak_extended_boundry import create_mp_slf, export_paths_to_shapefile
from waimak_extended_boundry import get_zone_array_index

def set_up_run_no_pumping_chch_selwyn(model_id, basedir, name):

    wells = get_base_well(model_id,True)
    wells = wells.loc[~((wells.zone == 's_wai') & (wells.use_type != 'injection'))]
    m = mod_gns_model(model_id,
                      name,
                      dir_path=os.path.join(basedir,name),
                      well={0: smt.convert_well_data_to_stresspd(wells)},
                      set_hdry=True)
    m.write_input()
    m.run_model()

    for zone in ['waimak','chch','selwyn']:
        print('starting zone {}'.format(zone))
        temp = np.array(smt.model_where(get_zone_array_index([zone])))
        particle_data = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(len(temp))
        particle_data['i0'][:] = temp[:,0]
        particle_data['j0'][:] = temp[:,1]
        particle_data['particlegroup'][:] = 1
        particle_data['groupname'][:] = 'part1'
        particle_data['label'][:] = 'test'
        mp = create_mp_slf(particle_data,m,os.path.join(basedir,name), mp_name=zone)
        mp.write_input()
        mp.write_name_file()
        mp.run_model()
        export_paths_to_shapefile(os.path.join(basedir,name, "{}.mppth".format(zone)),
                                  os.path.join(basedir,name, "{}_{}.shp".format(model_id, zone)))
if __name__ == '__main__':
    set_up_run_no_pumping_chch_selwyn('NsmcBase',"C:\Users\MattH\Desktop\simple_modpath_no_south_pumping",'simple_modpath_no_south_pumping')
