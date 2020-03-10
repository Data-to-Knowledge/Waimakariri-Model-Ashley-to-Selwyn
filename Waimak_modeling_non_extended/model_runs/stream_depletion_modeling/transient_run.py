"""
Author: matth
Date Created: 23/05/2017 12:25 PM
"""

from __future__ import division
import numpy as np
import users.MH.Waimak_modeling.model_tools as mt
import transient_inputs as tinputs
from copy import deepcopy
import flopy
import os
import sys
import pandas as pd



def setup_and_run_transient_multip(kwargs):
    """
    quick wrapper to make it easier to feed the below function to a multiprocessing script 
    :param kwargs: 
    :return: 
    """
    try:
        name, success = setup_and_run_transient(**kwargs)
        print(name, success)
        return name, success
    except Exception as val:
        print(val)
        sys.stdout.flush()

def setup_and_run_transient(name, base_dir, stress_vals, stress_to_month, wells_to_turn_off, ss, sy, silent=True,
                            solver = 'pcg', model='mf2005', start_heads=None,
                            model_version='m_strong_vert', well_data=None, rch_data=None, str_type='base'):
    """
    set up and run the stream depletion modflow models will write log to the base_dir
    :param name: name of the model
    :param base_dir: the directory to place the folder containing the model
    :param stress_vals: see stress_period_vals of mod_gns_model
    :param stress_to_month: a dictionary linking stress period value to the month e.g. {s_period: month, 0: 1} this must
                            be defined explicitly
    :param wells_to_turn_off: dictionary of list well names to turn off e.g. {s_period: list of wells}
                              if a stress period is absent the values from the previous stress period will be propogated
                              to the current stress period.  if want no wells to be turned off pass an empty list
    :param ss: specific storage shape (17,190,365) or float
    :param sy: specific yield shape (17,190,365) or float
    :param silent: passed to m.run_model.  if true Mirror modflow information to consol
    :param solver: default is pcg, only gmg, pcgn is also setup
    :param model: which version of modflow to use (default 2005, also supports nwt)
    :param model_version: which version of the waimak model to use
    :param well_data: a dataframe of well data to use in the simulation if none than uses the base well data
    :param rch_data: a array of recharge to use as the base values
    :param start_heads: the starting heads to use for the modflow model if None then defaults to the GNS provided starting heads
    :param str_type: 'base': use the GNS value, 'nat': as base, but without the WIL bywash into the cust.
    :return: 
    """
    # check inputs are dictionaries
    for input_arg in ['stress_vals', 'stress_to_month', 'wells_to_turn_off']:
        if not isinstance(eval(input_arg), dict):
            raise ValueError('incorrect input type for {} expected dict'.format(input_arg))

    nper = stress_vals['nper']

    # propigate last defined wells and stream depletion values to turn off forward
    if 0 not in wells_to_turn_off.keys():
        raise ValueError('no inital wells to turn of set.  if none set stress period 0 to an empty list')
    for per in range(nper):
        if per not in wells_to_turn_off.keys():
            wells_to_turn_off[per] = wells_to_turn_off[per-1]

    wells = {}
    rch = {}
    stream = {}
    if str_type == 'base':
        base_str = mt.get_base_str()
    elif str_type == 'nat':
        base_str = mt.get_base_str()
        base_str['flow'][827] = 0
    else:
        raise ValueError('unexpected argument for str_type: {}'.format(str_type))

    if rch_data is None:
        base_rch = mt.get_base_rch()
    elif isinstance(rch_data,np.ndarray):
        if rch_data.shape != (190,365):
            raise ValueError('incorrect shape for rch data: {}'.format(rch_data.shape))
        base_rch = rch_data

    if well_data is None:
        base_well = mt.get_original_well_data()
    elif isinstance(well_data,pd.DataFrame):
        if not {'row', 'col', 'layer', 'flux'}.issubset(set(well_data.keys())):
            raise ValueError('unexpected format for well_data missing keys: {}'.format({'row','col','layer','flux'}
                                                                                       -set(well_data.keys())))
        base_well = well_data

    else:
        raise ValueError('unexpected value for well_data: {}'.format(type(well_data)))

    rch_weighting = tinputs.get_mean_monthly_lsr()
    str_weighting = tinputs.get_monthly_str_scale()
    well_weighting = tinputs.get_mean_month_well_scaling()

    for sp in range(nper):

        month = stress_to_month[sp]  # get the month of the stress period
        if month is None:  # use the average value (to support steady state runs)
            stream[sp] = base_str
            rch[sp] = base_rch
            input_wells = deepcopy(base_well)
            if wells_to_turn_off[sp] != []:
                input_wells.flux[np.in1d(input_wells.index, wells_to_turn_off[sp])] = 0

            wells[sp] = mt.convert_well_data_to_stresspd(input_wells)
            continue

        # set up recharge
        rch[sp] = base_rch * rch_weighting.loc[month]

        changes = {  # set up stream
            'ashley': str_weighting['ashley'][month],
            # ashley northern tribs
            'glentui': str_weighting['a_tribs'][month],
            'garry': str_weighting['a_tribs'][month],
            'bullock': str_weighting['a_tribs'][month],
            'okuku': str_weighting['a_tribs'][month],
            'makerikeri': str_weighting['a_tribs'][month],
            'eyre': str_weighting['eyre'][month]}

        stream[sp] = mt.change_inflow(base_str, changes, set_or_mult='mult')

        # set up wells
        # apply seasonal weighting
        irr_idx = (base_well.use1 == 'irrigation') | (base_well.use2 == 'irrigation')
        input_wells = deepcopy(base_well)
        input_wells.flux[irr_idx] *= well_weighting['irrigation'].loc[month]
        input_wells.flux[
            ((base_well.use1 == 'public_water_supply') | (base_well.use2 == 'public_water_supply')) & (~irr_idx)] *= \
            well_weighting['public_water_supply'][month]

        #turn off selected wells
        input_wells.flux[np.in1d(input_wells.index, wells_to_turn_off[sp])] = 0
        wells[sp] = mt.convert_well_data_to_stresspd(input_wells)

    m = mt.wraps.mflow.mod_gns_model(name, '{}/{}'.format(base_dir, name),
                                     safe_mode=False,
                                     stress_period_vals=stress_vals,
                                     well=wells,
                                     drain=None,  # not modifying these stress period data
                                     recharge=rch,
                                     stream=stream,
                                     mt3d_link=False,
                                     start_heads=start_heads,
                                     model_version=model_version)

    # set specific storage and specific yield
    ss = np.atleast_1d(ss)
    sy = np.atleast_1d(sy)


    if len(ss) == 1:
        ss = np.ones((17,190,365)) * ss[0]

    if len(sy) == 1:
        sy = np.ones((17,190,365)) * sy[0]

    m.lpf.ss = flopy.utils.Util3d(m, m.lpf.ss.shape, m.lpf.ss.dtype, ss, m.lpf.ss.name)
    m.lpf.sy = flopy.utils.Util3d(m, m.lpf.sy.shape, m.lpf.sy.dtype, sy, m.lpf.sy.name)

    #set up re-wetting parameters and set any negative thickness to 1
    l1thick = m.dis.top.array - m.dis.botm.array[0]
    neg_idx = l1thick < 0
    top = m.dis.top.array
    top[neg_idx] += 50
    m.dis.top = flopy.utils.Util2d(m,m.dis.top.shape,m.dis.top.dtype,
                                    top,m.dis.top.name)

    m.lpf.laywet = flopy.utils.Util2d(m,m.lpf.laywet.shape,m.lpf.laywet.dtype,
                                      m.lpf.laytyp.array,m.lpf.laywet.name)
    m.lpf.ihdwet = 0
    m.lpf.iwetit = 1
    m.lpf.wetfct = 0.1
    m.lpf.wetdry = flopy.utils.Util3d(m,m.lpf.wetdry.shape,m.lpf.wetdry.dtype,
                                      np.abs(m.bas6.ibound.array),m.lpf.wetdry.name)

    if model.lower() == 'mf2005':
        pass
    elif model.lower() == 'mfnwt':
        mt.wraps.convert_to_nwt(m, iphdry=1)
        solver = 'nwt'

    else:
        raise ValueError('model {} not implemented'.format(model))
    # change set up for pCG
    m.pcg.mxiter = 200
    m.pcg.hclose = 1e-4
    m.pcg.rclose = 0.1
    if solver.lower() != 'pcg':
        m.remove_package('PCG')

    if solver.lower() == 'gmg':
        # change to GMG
        flopy.modflow.mfgmg.ModflowGmg(m,
                                       mxiter=300,
                                       iiter=1000, # might be able to be reduced
                                       iadamp=1,
                                       hclose=1e-04,
                                       rclose=0.1, #done to match pcg #1e-04,
                                       relax=1.0,
                                       ioutgmg=1,
                                       iunitmhc=None,
                                       ism=0,
                                       isc=1,
                                       damp=0.01,
                                       dup=0.4,
                                       dlow=0.01,
                                       chglimit=1.0,
                                       unitnumber=714)
    elif solver.lower() == 'pcgn':
        flopy.modflow.mfpcgn.ModflowPcgn(m,
                                         iter_mo=200, #lower?
                                         iter_mi=1000,
                                         close_r=1e-03,
                                         close_h=1e-03,
                                         relax=0.99,
                                         ifill=1,
                                         adamp=1,
                                         damp=1.0,
                                         damp_lb=0.001,
                                         rate_d=0.05,
                                         chglimit=0.0,
                                         acnvg=2,
                                         cnvg_lb=0.001,
                                         mcnvg=2,
                                         rate_c=-1.0,
                                         ipunit=0,
                                         unitnumber=714)
    elif solver.lower() == 'pcg':
        pass
    elif solver.lower() == 'nwt':
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
    else:
        raise ValueError('{} solver not implemented'.format(solver))
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
