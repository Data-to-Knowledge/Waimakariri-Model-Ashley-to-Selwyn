"""
Author: matth
Date Created: 11/05/2017 4:29 PM
"""

from __future__ import division
import users.MH.Waimak_modeling.model_tools as mt

# generate particles
wells = mt.get_original_well_data()
wells = wells[wells.flux < 0]
particle_data = mt.wraps.mpath.particle_loc_from_grid((wells['layer'],wells['row'],wells['col']))
particle_data = particle_data[particle_data['j0'] != 0]
particle_data = particle_data[particle_data['i0'] != 0]

# set up a modflow run
m = mt.wraps.mflow.import_gns_model('test',r"C:\Users\MattH\Desktop\test_dir")
m.write_input()
m.write_name_file()
m.run_model()


mp = mt.wraps.mpath.create_mp_slf(m, particle_data,
                                      prsity=0.15,
                                      prsityCB=0.15,
                                      mp_name=None,
                                      direction='backward',
                                      simulation_type='pathline',
                                      capt_weak_s=True,
                                      time_pts=1)
mp.write_input()
mp.write_name_file()
mp_success, buff1 = mp.run_model()

