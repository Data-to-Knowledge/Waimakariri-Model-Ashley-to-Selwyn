# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 5/12/2017 11:13 AM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.modpath_wrapper import export_paths_to_shapefile, create_mp_slf
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.base_modflow_wrapper import import_gns_model
import flopy
import numpy as np

if __name__ == '__main__':
    # run GNS model

    m = import_gns_model('StrOpt','modpath_base',r"C:\Users\MattH\Desktop\simple_modpath",False)
    m.write_name_file()
    m.upw.iphdry = 0  # hdry is -888.0
    m.write_input()
    m.run_model()

    # set up and run modpath models

    # waimak
    for zone in ['waimak','chch','selwyn']:
        print('starting zone {}'.format(zone))
        temp = np.array(smt.model_where(get_zone_array_index([zone])))
        particle_data = flopy.modpath.mpsim.StartingLocationsFile.get_empty_starting_locations_data(len(temp))
        particle_data['i0'][:] = temp[:,0]
        particle_data['j0'][:] = temp[:,1]
        particle_data['particlegroup'][:] = 1
        particle_data['groupname'][:] = 'part1'
        particle_data['label'][:] = 'test'
        mp = create_mp_slf(particle_data,m,r"C:\Users\MattH\Desktop\NsmcBase_simple_modpath\modpath_files",mp_name=zone)
        mp.write_input()
        mp.write_name_file()
        mp.run_model()
        export_paths_to_shapefile(r"C:\Users\MattH\Desktop\NsmcBase_simple_modpath\modpath_files\{}.mppth".format(zone),
                                  r"C:\Users\MattH\Desktop\NsmcBase_simple_modpath\{}.shp".format(zone))
