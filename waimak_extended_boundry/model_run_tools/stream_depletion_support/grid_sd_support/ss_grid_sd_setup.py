# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/10/2017 11:45 AM
"""

from __future__ import division
import env
import numpy as np
from copy import deepcopy
import flopy_mh as flopy
import os
import sys
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools import mod_gns_model, get_full_consent, get_race_data, zip_non_essential_files, modflow_converged
from traceback import format_exc
import pandas as pd
import pickle
from waimak_extended_boundry.model_run_tools.model_bc_data.drn_data import get_drn_no_ncarpet_spd



def setup_and_run_ss_grid_stream_dep_multip(kwargs):
    """
    quick wrapper to make it easier to feed the below function to a multiprocessing script
    :param kwargs:
    :return:
    """
    try:
        name, success = setup_and_run_ss_grid_stream_dep(**kwargs)
        print(name, success)
        return name, success
    except Exception as val:
        name = kwargs['name']
        success = format_exc().replace('\n', '')
    return name, success


def setup_and_run_ss_grid_stream_dep(model_id, name, base_dir, wells_to_turn_on,
                                     silent=True, start_heads=None):
    """
    set up and run the stream depletion modflow models will write log to the base_dir
    :param model_id: the model id for the NSMC realisation
    :param name: name of the model
    :param base_dir: the directory to place the folder containing the model
    :param wells_to_turn_on: dataframe of well  if none pass None
    :param silent: passed to m.run_model.  if true Mirror modflow information to consol
    :param start_heads: starting heads to pass to the model
    :return:
    """

    # check inputs are dictionaries
    for input_arg in ['wells_to_turn_on'] :
        if not isinstance(eval(input_arg), pd.DataFrame) and eval(input_arg) is not None:
            raise ValueError('incorrect input type for {} expected dataframe'.format(input_arg))

    wells = {}
    drns = get_drn_no_ncarpet_spd(model_id)
    base_well = get_race_data(model_id)
    full_consent = get_full_consent(model_id)
    # set up wells
    input_wells = deepcopy(base_well)
    if wells_to_turn_on is not None:
        input_wells = pd.concat((input_wells, wells_to_turn_on))

    wells[0] = smt.convert_well_data_to_stresspd(input_wells)

    m = mod_gns_model(model_id, name, '{}/{}'.format(base_dir, name),
                      safe_mode=False,
                      well=wells,
                      drain={0:drns},
                      mt3d_link=False,
                      start_heads=start_heads)

    # below included for easy manipulation
    flopy.modflow.mfnwt.ModflowNwt(m,
                                   headtol=1e-5,
                                   fluxtol=500,
                                   maxiterout=100,
                                   thickfact=1e-05,
                                   linmeth=1,
                                   iprnwt=0,
                                   ibotav=0,
                                   options='COMPLEX',
                                   Continue=False,
                                   dbdtheta=0.4,  # only when options is specified
                                   dbdkappa=1e-05,  # only when options is specified
                                   dbdgamma=0.0,  # only when options is specified
                                   momfact=0.1,  # only when options is specified
                                   backflag=1,  # only when options is specified
                                   maxbackiter=50,  # only when options is specified
                                   backtol=1.1,  # only when options is specified
                                   backreduce=0.7,  # only when options is specified
                                   maxitinner=50,  # only when options is specified
                                   ilumethod=2,  # only when options is specified
                                   levfill=5,  # only when options is specified
                                   stoptol=1e-10,  # only when options is specified
                                   msdr=15,  # only when options is specified
                                   iacl=2,  # only when options is specified
                                   norder=1,  # only when options is specified
                                   level=5,  # only when options is specified
                                   north=7,  # only when options is specified
                                   iredsys=0,  # only when options is specified
                                   rrctols=0.0,  # only when options is specified
                                   idroptol=1,  # only when options is specified
                                   epsrn=0.0001,  # only when options is specified
                                   hclosexmd=0.0001,  # only when options is specified
                                   mxiterxmd=50,  # only when options is specified
                                   unitnumber=714)

    # write inputs and run the model and write output to a log
    m.write_input()
    m.write_name_file()
    if silent:
        print('starting to run model {}'.format(name))
        sys.stdout.flush()
    success, buff = m.run_model(silent=silent, report=True)

    # write a log of the buffer
    log_dir = '{}/logging'.format(base_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log = '{}/{}_log.txt'.format(log_dir, name)
    buff = [e + '\n' for e in buff]
    with open(log, 'w') as f:
        f.writelines(buff)

    # get success and zip files I don't need for this analysis
    con = None
    if success:
        con = modflow_converged(os.path.join(m.model_ws, m.namefile.replace('.nam', '.list')))
        zip_non_essential_files(m.model_ws, include_list=False, other_files=['.sfo', '.ddn', '.hds'])
    if con is None:
        success = 'convergence unknown'
    elif con:
        success = 'converged'
    else:
        success = 'did not converge'
    return name, success


def grid_wells(flux, recalc=False):  # set up a grid
    pickle_path = os.path.join(smt.temp_pickle_dir, 'grid_sd_wells.p')

    if os.path.exists(pickle_path) and not recalc:
        outdata = pickle.load(open(pickle_path))
    else:
        base_dir = env.sci(
            "Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\grid_sd\inputs")
        layer_shapes = {0: 'layer1_2.shp',  # the file names are 1 indexed while they python layers are zero indexed
                        1: 'layer1_2.shp',
                        2: 'layer3_5_clip2.shp',
                        3: 'layer3_5_clip2.shp',
                        4: 'layer3_5_clip2.shp',
                        5: 'layer6_10_clip2.shp',
                        6: 'layer6_10_clip2.shp',
                        7: 'layer6_10_clip2.shp',
                        8: 'layer6_10_clip2.shp',
                        9: 'layer6_10_clip2.shp'}
        outdata = []
        for layer, path in layer_shapes.items():
            no_flow = smt.get_no_flow(layer)
            no_flow[no_flow < 0] = 0
            temp_array = smt.shape_file_to_model_array(os.path.join(base_dir, path), 'Id', True)
            temp_array[~no_flow.astype(bool)] = np.nan
            temp = np.array(smt.model_where(np.isfinite(temp_array)))

            rows, cols = temp[:, 0], temp[:, 1]
            layers = (np.zeros(rows.shape) + layer).astype(int)
            temp_df = pd.DataFrame({'row': rows, 'col': cols, 'layer': layers})
            outdata.append(temp_df)
        outdata = pd.concat(outdata).reset_index()
        elv_db = smt.calc_elv_db()
        for i in outdata.index:
            layer, row, col = outdata.loc[i, ['layer', 'row', 'col']].astype(int)
            x, y, z = smt.convert_matrix_to_coords(row, col, layer, elv_db)
            outdata.loc[i, 'mx'] = x
            outdata.loc[i, 'my'] = y
            outdata.loc[i, 'mz'] = z
            outdata.loc[i, 'depth'] = elv_db[0, row, col] - z
        pickle.dump(outdata, open(pickle_path, 'w'))

    outdata.loc[:, 'flux'] = flux
    outdata.loc[:, 'name'] = ['well_kijf_{}_{}_{}_{:.2f}'.format(k, i, j, f) for k, i, j, f in
                              outdata.loc[:, ['layer', 'row', 'col', 'flux']].itertuples(False, None)]
    outdata = outdata.set_index('name')
    return outdata


def get_base_grid_sd_path(model_id, recalc=False):
    path = os.path.join(smt.temp_pickle_dir, '{}_grid_sd_base'.format(model_id), '{}_grid_sd_base.cbc'.format(model_id))
    if os.path.exists(path) and not recalc:
        return path.replace('.cbc', '')

    name, success = setup_and_run_ss_grid_stream_dep(model_id=model_id,
                                     name='grid_sd_base',
                                     base_dir=smt.temp_pickle_dir,
                                     wells_to_turn_on=None)

    if not os.path.exists(path):
        raise ValueError('for some reason the cbc path did not write')

    if success != 'converged':
        os.remove(path) # to prevent it from returing the path on a future run
        raise ValueError('base model did not converge')

    return path.replace('.cbc', '')


if __name__ == '__main__':
    test = grid_wells(10, True)
    temp = smt.df_to_array(test, 'depth', True)
    for i in range(smt.layers):
        smt.plt_matrix(temp[i], title=i,base_map=True,alpha=1)
    import matplotlib.pyplot as plt

    plt.show()
    print('done')
