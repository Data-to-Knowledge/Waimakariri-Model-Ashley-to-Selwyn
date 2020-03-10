"""
Author: matth
Date Created: 4/05/2017 8:57 AM
"""

from __future__ import division
import flopy
import numpy as np
import pandas as pd
import time
import os
import pickle
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, sdp



def run_domestic_well_catch_model(m, name, weaksource, base_par_data):
    t = time.time()
    particle_data = base_par_data

    mp = mt.wraps.mpath.create_mp_slf(m, particle_data,
                                      prsity=0.15,
                                      prsityCB=0.15,
                                      mp_name=name,
                                      direction='backward',
                                      simulation_type='pathline',
                                      capt_weak_s=weaksource,
                                      time_pts=1)
    mp.write_input()
    mp.write_name_file()
    mp_success, buff1 = mp.run_model()

    if not mp_success:
        raise AssertionError('modpath model did not run')
    print('{} took {} minutes to run'.format(name, (time.time() - t) / 60))

def run_forward_simulation(m, name, weaksource):
    t = time.time()
    mp = mt.wraps.mpath.create_mp(m,
                                  (0, 0, 1),
                                  trackdir='forward',
                                  packages='RCH',
                                  simtype='pathline',
                                  time_pts=None,
                                  capt_weak_s=weaksource,
                                  prsity=0.15,
                                  prsityCB=0.15,
                                  ParticleColumnCount=4,
                                  ParticleRowCount=4,
                                  mp_name=name
                                  )
    mp.write_input()
    mp.write_name_file()
    mp_success, buff1 = mp.run_model()

    if not mp_success:
        raise AssertionError('modpath model did not run')

    print('{} took {} minutes to run'.format(name,(time.time()-t)/60))

def get_domestic_particles(recalc=False):
    pickle_path = '{}/inputs/pickeled_files/domestic_well_particles'.format(sdp)
    if os.path.exists(pickle_path) and not recalc:
        base_particle_data = pickle.load(open(pickle_path))
        return base_particle_data

    uses = ['domestic_supply', 'domestic_stockwater']
    well_data = mt.get_all_well_data()
    well_data_in = well_data[(np.in1d(well_data['use1'], uses)) | (np.in1d(well_data['use2'], uses))]

    # put particle rings (5 particles) around each well at the top/middle/bottom of the screens 25m away
    base_particle_data = mt.wraps.mpath.particle_ring_around_screens(25, 25,
                                                                     xs=well_data_in['x'], ys=well_data_in['y'],
                                                                     screen_zs=well_data_in['screen_z'],
                                                                     labels=well_data_in.index.str.replace('/','_'),
                                                                     forgive=True)
    # remove wells/particles which are on the boundry of the model(sides) as this causes problems in the endpoint file
    base_particle_data = base_particle_data[base_particle_data['j0'] != 0]
    base_particle_data = base_particle_data[base_particle_data['i0'] != 0]

    pickle.dump(base_particle_data, open(pickle_path,'w'))
    return base_particle_data

def run_models(raise_except=True):
    model_path = '{}/well_flow_paths'.format(base_mod_dir)
    m = mt.wraps.mflow.import_gns_model('flow_paths', model_path)
    m.write_input()
    m.write_name_file()
    m_success, buff = m.run_model()

    if not m_success:
        raise AssertionError('modflow model did not run')

    print('get domestic particle data')
    base_particle_data = get_domestic_particles()
    print('running domestic simulations')
    if raise_except:
        run_domestic_well_catch_model(m, 'domestic_wells_weak_s', True, base_particle_data)
        run_domestic_well_catch_model(m, 'domestic_wells_Strong_s', False, base_particle_data)
    else:
        try:
            run_domestic_well_catch_model(m, 'domestic_wells_weak_s', True, base_particle_data)
        except Exception as val:
            print('domestic_wells_weak_s ', val)

        try:
            run_domestic_well_catch_model(m, 'domestic_wells_Strong_s', False, base_particle_data)
        except Exception as val:
            print('domestic_wells_Strong_s', val)

    # run water supply simulations

    wdc_path = "{}/inputs/wells/WDC groundwater takes.xlsx".format(sdp)

    temp_data = pd.read_excel(wdc_path, sheetname='Final Data')
    well_data = mt.get_all_well_data()
    wdc_idx = np.in1d(well_data.index, np.array(temp_data['Well_no']))
    wdc_wells = well_data[wdc_idx]

    print('running forward WDC cells')
    # run a forward simulations by placing a 16 particles in each recharge cell
    if raise_except:
        run_forward_simulation(m, 'wdc_forward_wells_weak_s', True)
        run_forward_simulation(m, 'wdc_forward_wells_strong_s', False)
    else:
        try:
            run_forward_simulation(m, 'wdc_forward_wells_weak_s', True)
        except Exception as val:
            print('wdc_forward_wells_weak_s', val)

        try:
            run_forward_simulation(m, 'wdc_forward_wells_strong_s', False)
        except Exception as val:
            print('wdc_forward_wells_strong_s', val)

    print('getting WDC particle data')
    # run backwards simulations
    # put particle rings (25 particles) around each well at the top/middle/bottom of the screens 25m away
    wdc_particle_data = mt.wraps.mpath.particle_ring_around_screens(25, 25, wdc_wells['x'], wdc_wells['y'],
                                                                    wdc_wells['screen_z'],
                                                                    labels=wdc_wells.index, forgive=True)

    print ('running backward WDC models')
    if raise_except:
        run_domestic_well_catch_model(m, 'wdc_back_wells_weak_s', True, wdc_particle_data)
        run_domestic_well_catch_model(m, 'wdc_back_wells_strong_s', False, wdc_particle_data)
    else:
        try:
            run_domestic_well_catch_model(m, 'wdc_back_wells_weak_s', True, wdc_particle_data)
        except Exception as val:
            print('wdc_back_wells_weak_s', val)

        try:
            run_domestic_well_catch_model(m, 'wdc_back_wells_strong_s', False, wdc_particle_data)
        except Exception as val:
            print('wdc_back_wells_strong_s', val)

    print('finished running modpath models')


def gen_backward_output(path):
    # load
    epo = flopy.utils.EndpointFile(path)
    part_ends = epo.get_alldata()

    # need to figure out shapes should be number of particles
    array_shp = (190,365)
    # initalize array
    outdata = np.zeros(array_shp)

    # get set of wells from labels
    wells = set(part_ends['label'])

    for well in wells:
        rows = part_ends['i'][part_ends['label'] == well]
        cols = part_ends['j'][part_ends['label'] == well]

        outdata[rows,cols] += 1 # check this thourouly

    return outdata

def gen_forward_output(path): # the labels are not working only keeping the first bit replace / with _ and see if that sorts it
    well_data = mt.get_all_well_data()
    wdc_path = "{}/inputs/wells/WDC groundwater takes.xlsx".format(sdp)

    temp_data = pd.read_excel(wdc_path, sheetname='Final Data')

    wdc_idx = np.in1d(well_data.index, np.array(temp_data['Well_no']))
    wdc_wells = well_data[wdc_idx]
    epo = flopy.utils.EndpointFile(path)

    array_shp = (190, 365)
    outdata = np.zeros(array_shp)
    for well in wdc_wells.index:
        layer = wdc_wells.loc[well, 'layer']
        row = wdc_wells.loc[well, 'row']
        col = wdc_wells.loc[well, 'col']
        rch_start = epo.get_destination_endpoint_data(dest_cells=[(layer,row,col)])
        rows = rch_start['i']
        cols = rch_start['j']
        outdata[rows, cols] += 1  # check this thourouly

    return outdata

def get_all_mod_output():
    model_path = '{}/well_flow_paths'.format(base_mod_dir)
    back_models = ['domestic_wells_weak_s', 'domestic_wells_Strong_s',
                   'wdc_back_wells_weak_s', 'wdc_back_wells_strong_s']
    sources = {}
    for model in back_models:
        mpath = '{}/{}.mpend'.format(model_path, model)
        sources[model] = gen_backward_output(mpath)

    forward_mod = ['wdc_forward_wells_weak_s', 'wdc_forward_wells_strong_s']
    for model in forward_mod:
        mpath = '{}/{}.mpend'.format(model_path, model)
        sources[model] = gen_forward_output(mpath)

    return sources


def main(recalc = False):
    run_models(raise_except=False)

    outdir = '{}/well_flow_paths/output'.format(base_mod_dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    output = get_all_mod_output()

    for key in output.keys():
        fig, ax = mt.plt_matrix(output[key])
        fig.set_size_inches(18.5, 10.5)
        fig.savefig('{}/{}.png'.format(outdir,key))

if __name__ == '__main__':
    #main() # fix after debug

    outdir = '{}/well_flow_paths/output'.format(base_mod_dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    output = get_all_mod_output()

    for key in output.keys():
        fig, ax = mt.plt_matrix(output[key])
        fig.set_size_inches(18.5, 10.5)
        fig.savefig('{}/{}.png'.format(outdir, key))