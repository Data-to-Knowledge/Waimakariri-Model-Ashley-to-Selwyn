"""
Author: matth
Date Created: 18/04/2017 8:11 AM
Date Modified: 7/09/2017 2:19 PM
"""

from __future__ import division
from future.builtins import input
from core import env
import numpy as np
import flopy
import os
import shutil
import warnings
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from copy import deepcopy
from users.MH.Waimak_modeling.supporting_data_path import sdp

org_data_dir = "{}/from_GNS".format(sdp)


# todo as part of this make another function which gets the name file (or copies the files across)

def get_model(model_id):
    """
    load a flopy model instance of the model id
    :param model_id:
    :return: m flopy model instance
    """
    # todo check well packaged loads appropriately! yep this is a problem because we define the iface value,
    # todo but it isn't used, remove iface manually?
    # todo if it doesn't exist I will need to add the options flags to SFR package manually?
    m = None
    raise NotImplementedError()  # todo


def import_gns_model(model_id, name, dir_path, safe_mode=True, mt3d_link=False, exe_path=None):
    """
    sets up a model object for a steady state run with the base modflow files from GNS. This serves as a base to be
    modified for other  model runs. DOES NOT Write input Files, name file, or RUN MODFLOW.
    :param model_id: the id linking to the realisation
    :param name: the model name
    :param dir_path: path to the directory to create modflow files it does not need to exist
    :param safe_mode: if true ask for confirmation before deleting dirpath
    :param model_version: the model to use see supporting data paths for avalible versions. default to m_strong_vert
                          for backwards compatability
    :param exe_path: path to the modflow NWT exicutable
    :return: model object
    """
    # set name to incorporate the model version as the start of the name
    name = '{}-{}'.format(model_id, name)
    dir_path = '{}/{}-{}'.format(os.path.dirname(dir_path), model_id, os.path.basename(dir_path))
    # remove all previous files in the directory
    if os.path.exists(dir_path):
        if safe_mode:
            cont = input(
                'create_base_modflow_files will delete the directory:\n {} \n continue y/n\n'.format(dir_path)).lower()
            if cont == 'n':
                raise ValueError('script aborted so as not to overwrite {}'.format(dir_path))

        shutil.rmtree(dir_path)
        os.makedirs(dir_path)
    else:
        os.makedirs(dir_path)

    # load modflow and check certain packages exist only load certain packages so that others can be loaded explicitly
    # m = get_model(model_id) #todo reinstate later
    m = flopy.modflow.Modflow.load(
        r"C:\Users\MattH\Desktop\to_test_write\m_ex_bd_va-with_n_carpet\m_ex_bd_va-to_test_load.nam")  # todo only a test

    # change the naming and paths to files (essentially duplicate all data)
    m._set_name(name)
    m.change_model_ws(dir_path)
    if exe_path is None:
        m.exe_name = "{}/models_exes/MODFLOW-NWT_1.1.2/MODFLOW-NWT_1.1.2/bin/MODFLOW-NWT_64.exe".format(sdp)
    else:
        m.exe_name = exe_path

    # remove LMT package (will recreate if needed)
    if 'LMT6' in m.get_package_list():
        m.remove_package('LMT6')

    # if needed create MT3D link files
    if mt3d_link:
        mt3d_outunit = 54
        mt3d_outname = '{}_mt3d_link.ftl'.format(m.name)
        link = flopy.modflow.ModflowLmt(m, output_file_name=mt3d_outname, output_file_unit=mt3d_outunit, unitnumber=21)
        m.add_output(mt3d_outname, mt3d_outunit, True, 'LMT')

    return m


def mod_gns_model(model_id, name, dir_path, safe_mode=True, stress_period_vals=None, well=None, drain=None,
                  recharge=None,
                  stream=None, mt3d_link=False, start_heads=None, exe_path=None):
    """
    modify the gns model by changing the stress period data for one of the following packages WEL, DRN, RCH, STR, DIS

    :param model_id: the link to the model realisations
    :param name: name of the model
    :param dir_path: the working directory of the model if it doesn't exist it's created else it's deleted and created
    :param safe_mode: if True it requires user input to delete dir_path
    :param stress_period_vals: default None or dictionary these are the values to set the model stress period to if 
                               None then the default stress period values from GNS are used [1,1.,1,True,1.]
                               dictionary contains: 
                                 nper: number of stress periods (must be integer)
                                 perlen: the length of each stress period (as value, list, or array of floats)
                                 nstp: the number of steps in each stress period ( as value, list or array of ints)
                                 steady: stress period steady state (bool or list array of bool)
                                 oc_stress_per_data: None or oc stress period data (see flopy documentation for style)
                                                     if None or not present all period and steps will handled by
                                                     ['save head', 'save drawdown', 'save budget', 'print budget']
                                 tsmult: the timestep multiplier (float or list or array of floats) 
    :param well: dictionary format stress period data for the WEL package (see flopy.ModflowWell) or None and use
           default GNS data
    :param drain: dictionary format stress period data for the DRN package (see flopy.ModflowWell) or None and use
           default GNS data
    :param recharge: dictionary format stress period data for the RCH package (see flopy.ModflowWell) or None and use
           default GNS data
    :param stream: tuple or list of (segment_data, reach_data) for the SFR package or None and use
           default model data.  Varying stress period data is not yet implemented.
    :param mt3d_link: boolean if true write a MT3D-link module
    :param start_heads: None or array for the starting heads of the simulation.
    :param exe_path: path to the modflow NWT exicutable
    :return: model
    """

    # check inputs
    if stress_period_vals is not None:
        stress_period_vals = _check_stress_period_values(stress_period_vals)

    # check stress period data and raise warnings if it is not the same length as the number of periods
    if stress_period_vals is None:
        nper = 1
    else:
        nper = stress_period_vals['nper']
    print_lg_warning = False
    if nper == 1:
        pass  # do not raise stress period data warnings for 1 stress period models
    else:
        for var_name in ['well', 'drain', 'recharge', 'stream']:
            if eval(var_name) is None:
                print_lg_warning = True
                warnings.warn('using default GNS values for {}. '
                              'This only has data for one stress period'.format(var_name))
            elif len(eval(var_name)) != nper:
                print_lg_warning = True
                warnings.warn('{} has stress period data for {} out '
                              'of {} stress periods'.format(var_name, len(eval(var_name)), nper))
    if print_lg_warning:
        warnings.warn('One or more stress period data warnings: \n'
                      'Note that if the number of lists is smaller than the number of stress periods,\n'
                      'then the last list of wells will apply until the end of the simulation.\n'
                      'Full details of all options to specify stress_period_data can be found in the\n'
                      'flopy3 boundaries Notebook in the basic subdirectory of the examples directory')

    for var_name in ['well', 'drain', 'recharge', 'stream']:
        if eval(var_name) is None:
            pass
        elif len(eval(var_name)) > nper:
            raise ValueError('{} has more stress periods supplied than present in model'.format(var_name))

    # import model and change stress period if needed
    m = import_gns_model(model_id, name, dir_path, safe_mode=safe_mode, mt3d_link=mt3d_link, exe_path=exe_path)

    if stress_period_vals is not None:
        change_stress_period_settings(m, stress_period_vals)

    # add well stress period data
    if well is not None:
        wel = flopy.modflow.ModflowWel(m, ipakcb=740, stress_period_data=well, unitnumber=709)

    # add drain stress period data
    if drain is not None:
        drn = flopy.modflow.ModflowDrn(m,
                                       ipakcb=740,
                                       stress_period_data=drain,
                                       unitnumber=710,
                                       options=[])

    # add recharge stress period data
    if recharge is not None:
        rch = flopy.modflow.ModflowRch(m, nrchop=3, ipakcb=740, rech=recharge, unitnumber=716)

    # add stream stress period data
    if stream is not None:
        if (not isinstance(stream, tuple) and not isinstance(stream, list)) and len(stream) != 2:
            raise ValueError('incorrect input type for stream')
        segment_data = stream[0]
        reach_data = stream[1]
        if isinstance(segment_data, dict) or isinstance(reach_data, dict):
            raise NotImplementedError('varying SFR stress period data not yet implemented')
        sfr = flopy.modflow.ModflowSfr2(m,
                                        nstrm=len(reach_data),
                                        nss=len(segment_data),
                                        nsfrpar=0,  # this is not implemented in flopy not using
                                        nparseg=0,  # this is not implemented in flopy not using
                                        const=86400,
                                        dleak=0.0001,  # default, likely ok
                                        ipakcb=740,
                                        istcb2=0,
                                        isfropt=1,  # Not using as nstrm >0 therfore no Unsaturated zone flow
                                        nstrail=10,  # not using as no unsaturated zone flow
                                        isuzn=1,  # not using as no unsaturated zone flow
                                        nsfrsets=30,  # not using as no unsaturated zone flow
                                        irtflg=0,  # no transient sf routing
                                        numtim=2,  # Not using transient SFR
                                        weight=0.75,  # not using if irtflg = 0 (no transient SFR)
                                        flwtol=0.0001,  # not using if irtflg = 0 (no transient SFR)
                                        reach_data=reach_data,
                                        segment_data=segment_data,
                                        channel_geometry_data=None,  # using rectangular profile so not used
                                        channel_flow_data=None,  # icalc = 1 or 2 so not using
                                        dataset_5={0: [len(segment_data), 0, 0]},
                                        reachinput=True,  # using for reach input
                                        transroute=False,  # no transient sf routing
                                        tabfiles=False,  # not using
                                        tabfiles_dict=None,  # not using
                                        unit_number=717)

        raise NotImplementedError()

    if start_heads is not None:
        temp = m.get_package('BAS6').strt
        m.get_package('BAS6').strt = flopy.utils.Util3d(m, temp.shape, temp.dtype, start_heads, temp.name)
    return m


def _check_stress_period_values(spv):
    """
    check that the provided parameters are in the correct formats for a description of parameter see DIS object of flopy
    :param stress_period_values: dictionary of: 
                                  nper: number of stress periods (must be integer)
                                  perlen: the length of each stress period (as value, list, or array of floats)
                                  nstp: the number of steps in each stress period ( as value, list or array of ints)
                                  steady: stress period steady state (bool or list array of bool)
                                  oc_stress_per_data: None of oc stress period data (see flopy documentation for style)
                                                     if None or not persent all period and steps will handled by
                                                     ['save head', 'save drawdown', 'save budget', 'print budget']
                                  tsmult: the timestep multiplier (float or list or array of floats) 
    :return: stress_period_values formatted properly
    """
    # format most input as arrays
    spv['tsmult'] = np.atleast_1d(spv['tsmult'])
    spv['perlen'] = np.atleast_1d(spv['perlen'])
    spv['nstp'] = np.atleast_1d(spv['nstp'])
    spv['steady'] = np.atleast_1d(spv['steady'])

    # check input variable types
    if not isinstance(spv['nper'], int):
        raise ValueError('nper must be integer')

    if spv['steady'].dtype != bool:
        try:
            spv['steady'] = spv['steady'].astype(bool)
        except:
            raise ValueError('steady must be boolean or transmutable to boolean')

    try:
        spv['perlen'] = spv['perlen'].astype(float)
    except:
        raise ValueError('perlen must be castable to float')
    try:
        spv['tsmult '] = spv['tsmult'].astype(float)
    except:
        raise ValueError('tsmult must be castable to float')

    try:
        if all(spv['nstp'].astype(int) / spv['nstp'].astype(float) == 1):
            spv['nstp'] = spv['nstp'].astype(int)
    except:
        raise ValueError('nstp must be castable to int')

    if not all(spv['nstp'].astype(int) / spv['nstp'].astype(float) == 1):
        raise ValueError('expected whole number for nstp')

    # check lengths and propigate single values

    for name in ['perlen', 'nstp', 'steady', 'tsmult']:
        if spv[name].ndim != 1:
            raise ValueError('{} must be 1-d')
        elif len(spv[name]) != spv['nper']:
            if len(spv[name]) == 1:
                spv[name] = (np.ones((spv['nper'],)) * spv[name]).astype(spv[name].dtype)
            else:
                raise ValueError('{} must have the same length as the number of periods ({})'.format(name, spv['nper']))

    return spv


def change_stress_period_settings(m, spv):
    """
    changes a model dis object to a new set of stress period values and run a check
    For more info on variables see flopy.modflow.ModflowDis
    :param m: modflow model
    :param  spv: dictionary of: nper: number of stress periods (must be integer)
                                perlen: the length of each stress period (as value, list, or array of floats)
                                nstp: the number of steps in each stress period ( as value, list or array of ints)
                                steady: stress period steady state (bool or list array of bool)
                                oc_stress_per_data: None of oc stress period data (see flopy documentation for style)
                                                    if None or not persent all period and steps will handled by
                                                    ['save head', 'save drawdown', 'save budget', 'print budget']
                                tsmult: the timestep multiplier (float or list or array of floats) 
    :return:
    """
    # takes the model object (as input) and changes all of the stress period conditions
    # does not change stress period data e.g. well stress peiod data
    spv = _check_stress_period_values(spv)
    dis = m.get_package('dis')
    nper = spv['nper']
    if 'oc_stress_per_data' not in spv.keys():
        oc_stress_per_data = None
    else:
        oc_stress_per_data = spv['oc_stress_per_data']
    dis.nper = nper
    dis.perlen = flopy.utils.Util2d(m, (nper,), np.float32, spv['perlen'], 'perlen')
    dis.nstp = flopy.utils.Util2d(m, (nper,), int, spv['nstp'], 'nstp')
    dis.steady = flopy.utils.Util2d(m, (nper,), bool, spv['steady'], 'steady')
    dis.tsmult = flopy.utils.Util2d(m, (nper,), np.float32, spv['tsmult'], 'tsmult')
    dis.check()
    if oc_stress_per_data is None:
        m.oc.stress_period_data = {(0, 0): ['save head', 'save drawdown', 'save budget', 'print budget']}
    else:
        m.oc.stress_period_data = oc_stress_per_data


if __name__ == '__main__':
    import_gns_model(1, 'test', r"C:\Users\MattH\Desktop\test", False, True)
    print('done')  # todo check modGNS model (though it should be pretty good)
