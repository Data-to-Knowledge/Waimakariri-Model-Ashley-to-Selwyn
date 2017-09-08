"""
Author: matth
Date Created: 7/09/2017 3:43 PM
"""

from __future__ import division
import numpy as np
from copy import deepcopy
import flopy
import os
import sys
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools import mod_gns_model,get_max_rate, get_full_consent, get_race_data


def setup_and_run_stream_dep_multip(kwargs):
    """
    quick wrapper to make it easier to feed the below function to a multiprocessing script
    :param kwargs:
    :return:
    """
    try:
        name, success = setup_and_run_stream_dep(**kwargs)
        print(name, success)
        return name, success
    except Exception as val:
        with open('{}/error_log.txt'.format(kwargs['base_dir']),mode='a') as f:
            f.write('{}: {}\n'.format(kwargs['name'],val))
        print(val)
        sys.stdout.flush()

def setup_and_run_stream_dep(model_id, name, base_dir, stress_vals, wells_to_turn_on, ss, sy,
                             silent=True, start_heads=None, sd_7_150='sd150'):
    """
    set up and run the stream depletion modflow models will write log to the base_dir
    :param model_id: the model id for the NSMC realisation
    :param name: name of the model
    :param base_dir: the directory to place the folder containing the model
    :param stress_vals: see stress_period_vals of mod_gns_model
    :param wells_to_turn_on: dictionary of list well names to turn on e.g. {s_period: list of wells}
                              if a stress period is absent the values from the previous stress period will be propogated
                              to the current stress period.  if want no wells to be turned off pass an empty list
    :param ss: specific storage shape (k,i,j) or float
    :param sy: specific yield shape (k,i,j) or float
    :param silent: passed to m.run_model.  if true Mirror modflow information to consol
    :param sd_7_150: sd150: annual take over 5 stress periods
                     sd7: use the max rate if no max rate set as annual take over 150 days
    :return:
    """
    # check inputs are dictionaries
    for input_arg in ['stress_vals', 'stress_to_month', 'wells_to_turn_on']:
        if not isinstance(eval(input_arg), dict):
            raise ValueError('incorrect input type for {} expected dict'.format(input_arg))

    nper = stress_vals['nper']

    # propigate last defined wells and stream depletion values to turn off forward
    if 0 not in wells_to_turn_on.keys():
        raise ValueError('no inital wells to turn of set.  if none set stress period 0 to an empty list')
    for per in range(nper):
        if per not in wells_to_turn_on.keys():
            wells_to_turn_on[per] = wells_to_turn_on[per - 1]

    wells = {}
    rch = {}
    stream = {}

    base_well = get_race_data()
    if sd_7_150 == 'sd150':
        full_consent = get_full_consent()
    elif sd_7_150 == 'sd7':
        sd7_flux = get_max_rate()
    for sp in range(nper):
        # set up recharge
        rch[sp] = 0

        # set up wells
        input_wells = deepcopy(base_well)
        for well in wells_to_turn_on[sp]:
            if sd_7_150 == 'sd150':
                add_well = full_consent.loc[well] #todo this needs to be handled
                add_well.loc['flux'] *= 12/5 #todo EAV check with this
                input_wells.loc[well] = add_well
            elif sd_7_150 == 'sd7':
                add_well = sd7_flux.loc[well]
                input_wells.loc[well] = add_well
            else:
                raise ValueError('unexpected argument for sd_7_150: {}'.format(sd_7_150))

        wells[sp] = smt.convert_well_data_to_stresspd(input_wells)

    m = mod_gns_model(model_id, name, '{}/{}'.format(base_dir, name),
                                     safe_mode=False,
                                     stress_period_vals=stress_vals,
                                     well=wells,
                                     drain=None,  # not modifying these stress period data
                                     recharge=rch,
                                     stream=None,
                                     mt3d_link=False,
                                     start_heads=start_heads)

    # set specific storage and specific yield
    ss = np.atleast_1d(ss)
    sy = np.atleast_1d(sy)


    if len(ss) == 1:
        ss = np.ones((smt.layers,smt.rows,smt.cols)) * ss[0]

    if len(sy) == 1:
        sy = np.ones((smt.layers,smt.rows,smt.cols)) * sy[0]

    m.upw.ss = flopy.utils.Util3d(m, m.lpf.ss.shape, m.lpf.ss.dtype, ss, m.lpf.ss.name)
    m.upw.sy = flopy.utils.Util3d(m, m.lpf.sy.shape, m.lpf.sy.dtype, sy, m.lpf.sy.name)

    #todo do I need to include this?
    flopy.modflow.mfnwt.ModflowNwt(m,
                                   headtol=0.01,
                                   fluxtol=500,
                                   maxiterout=100,
                                   thickfact=1e-05,
                                   linmeth=1,
                                   iprnwt=0,
                                   ibotav=0,
                                   options='COMPLEX',
                                   Continue=False,
                                   dbdtheta=0.4,
                                   dbdkappa=1e-05,
                                   dbdgamma=0.0,
                                   momfact=0.1,
                                   backflag=1,
                                   maxbackiter=50,
                                   backtol=1.1,
                                   backreduce=0.7,
                                   maxitinner=50,
                                   ilumethod=2,
                                   levfill=5,
                                   stoptol=1e-10,
                                   msdr=15,
                                   iacl=2,
                                   norder=1,
                                   level=5,
                                   north=7,
                                   iredsys=0,
                                   rrctols=0.0,
                                   idroptol=1,
                                   epsrn=0.0001,
                                   hclosexmd=0.0001,
                                   mxiterxmd=50,
                                   unitnumber=714)
    # write inputs and run the model and write output to a log
    m.write_input()
    m.write_name_file()
    if silent:
        print('starting to run model {}'.format(name))
        sys.stdout.flush()
    success, buff = m.run_model(silent=silent, report=True)
    log_dir = '{}/logging'.format(base_dir)
    if not os.path.exists (log_dir):
        os.makedirs(log_dir)
    log = '{}/{}_log.txt'.format(log_dir, name)
    buff = [e + '\n' for e in buff]
    with open(log, 'w') as f:
        f.writelines(buff)
    return name, success
#todo this needs debugging