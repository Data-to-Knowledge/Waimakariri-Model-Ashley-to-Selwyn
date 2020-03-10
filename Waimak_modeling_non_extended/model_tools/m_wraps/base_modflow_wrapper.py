"""
Author: matth
Date Created: 18/04/2017 8:11 AM
"""

from __future__ import division
from future.builtins import input
import numpy as np
import flopy_mh as flopy
import os
import shutil
import warnings
from copy import deepcopy
from users.MH.Waimak_modeling.supporting_data_path import sdp, get_org_mod_path
from users.MH.Waimak_modeling.model_tools.well_values import get_original_well_data, convert_well_data_to_stresspd

org_data_dir = "{}/from_GNS".format(sdp)
def import_gns_model(name, dir_path, safe_mode=True, mt3d_link=False, model_version='m_strong_vert'):
    """
    sets up a model object for a steady state run with the base modflow files from GNS. This serves as a base to be
    modified for other  model runs. DOES NOT Write input Files, name file, or RUN MODFLOW.
    :param name: the model name
    :param dir_path: path to the directory to create modflow files it does not need to exist
    :param safe_mode: if true ask for confirmation before deleting the directory
    :param model_version: the model to use see supporting data paths for avalible versions. default to m_strong_vert
                          for backwards compatability
    :return: model object
    """
    # set name to incorporate the model version as the start of the name
    name = '{}-{}'.format(model_version, name)
    dir_path = '{}/{}-{}'.format(os.path.dirname(dir_path), model_version, os.path.basename(dir_path))
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
    org_mod_path = get_org_mod_path(model_version)
    # load modflow and check certain packages exist only load certain packages so that others can be loaded explicitly
    m = flopy.modflow.Modflow.load('{}.nam'.format(org_mod_path), model_ws=os.path.dirname(org_mod_path), forgive=False)

    print('the model has been loaded with the following packages: {}'.format(m.get_package_list()))
    # change the naming and paths to files (essentially duplicate all data)
    m._set_name(name)
    m.change_model_ws(dir_path)
    m.exe_name = "{}/models_exes/mf2005.exe".format(sdp)

    # remove all old output files
    units = deepcopy(m.output_units)
    for u in units:
        m.remove_output(unit=u)

    # set package output units for those that need to be specified
    m.get_package('DRN').ipakcb = 741
    m.get_package('RCH').ipakcb = 742
    m.get_package('STR').ipakcb = 743 # between reach and model cells
    m.get_package('STR').istcb2 = 744 # between reach and reach

    # remove LMT package (will recreate if needed)
    if 'LMT6' in m.get_package_list():
        m.remove_package('LMT6')


    # add output to name file
    fnames = [m.name + e for e in ['.hds', '.ddn', '.cbb', '.cbd', '.crc', '.cstm', '.csts', '.cbw']]  # output extension
    funits = [30, 31, 740, 741, 742, 743, 744, 745] #fortran unit
    fbinflag = [True, True, True, True, True, True, True, True, True] #is binary
    fpackage = [[], [], ['LPF'], ['DRN'], ['RCH'], ['STR'], ['STR'], ['WEL']]

    # if needed create MT3D link files
    if mt3d_link:
        mt3d_outunit = 54
        mt3d_outname = '{}_mt3d_link.ftl'.format(m.name)
        link = flopy.modflow.ModflowLmt(m, output_file_name=mt3d_outname, output_file_unit=mt3d_outunit,unitnumber=21)

        fnames.append(mt3d_outname)
        funits.append(mt3d_outunit)
        fbinflag.append(True)
        fpackage.append(['LMT'])

    for fn, fu, fb, fp, in zip(fnames, funits, fbinflag, fpackage):
        m.add_output(fn,fu,fb,fp)

    print('loading well package')
    # explicitly add WEL package as it is acting up on the original load
    # load stress period data
    well_data = get_original_well_data()
    spd = convert_well_data_to_stresspd(well_data)
    wel = flopy.modflow.ModflowWel(m,ipakcb=745,stress_period_data = {0:spd},options=['AUX IFACE'])

    return m

def mod_gns_model(name, dir_path, safe_mode=True, stress_period_vals=None, well=None, drain=None, recharge=None,
                  stream=None, mt3d_link=False, start_heads=None, model_version='m_strong_vert'):
    """
    modify the gns model by changing the stress period data for one of the following packages WEL, DRN, RCH, STR, DIS

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
    :param stream: dictionary format stress period data for the STR package (see flopy.ModflowWell) or None and use
           default GNS data
    :param mt3d_link: boolean if true write a MT3D-link module
    :param start_heads: None or array for the starting heads of the simulation.
    :param model_version: the model to use see supporting data paths for avalible versions. default to m_strong_vert
                          for backwards compatability
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
        pass # do not raise stress period data warnings for 1 stress period models
    else:
        for var_name in ['well', 'drain', 'recharge', 'stream']:
            if eval(var_name) is None:
                print_lg_warning = True
                warnings.warn('using default GNS values for {}. '
                              'This only has data for one stress period'.format(var_name))
            elif len(eval(var_name)) != nper:
                print_lg_warning = True
                warnings.warn('{} has stress period data for {} out '
                              'of {} stress periods'.format(var_name,len(eval(var_name)),nper))
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
    m = import_gns_model(name, dir_path, safe_mode=safe_mode, mt3d_link=mt3d_link, model_version=model_version)

    if stress_period_vals is not None:
        change_stress_period_settings(m, stress_period_vals)

    # add well stress period data
    if well is not None:
        wel = flopy.modflow.ModflowWel(m, ipakcb=745, stress_period_data=well, options=['AUX IFACE'])

    # add drain stress period data
    if drain is not None:
        drn = flopy.modflow.ModflowDrn(m,
                                       ipakcb=741,
                                       stress_period_data=drain,
                                       unitnumber=710,
                                       options=[])

    # add recharge stress period data
    if recharge is not None:
        rch = flopy.modflow.ModflowRch(m, nrchop=3, ipakcb=742, rech=recharge, unitnumber=716)


    # add stream stress period data
    if stream is not None:
        seg = m.str.segment_data
        for key in stream.keys():
            if key == 0:
                continue
            seg[key] = seg[0]
        str_ = flopy.modflow.mfstr.ModflowStr(m,
                                               mxacts=1217,
                                               nss=44,
                                               ntrib=3,
                                               ndiv=0,
                                               icalc=1,
                                               const=86400.0,
                                               ipakcb=743,
                                               istcb2=744,
                                               dtype=None,
                                               stress_period_data=stream,
                                               segment_data=seg,
                                               extension='str',
                                               unitnumber=717)

    if start_heads is not None:
        temp = m.get_package('BAS6').strt
        m.get_package('BAS6').strt = flopy.utils.Util3d(m,temp.shape,temp.dtype,start_heads, temp.name)
    return m

def _check_stress_period_values (spv):
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
        spv['tsmult ']= spv['tsmult'].astype(float)
    except:
        raise ValueError('tsmult must be castable to float')

    try:
        if all(spv['nstp'].astype(int)/spv['nstp'].astype(float) == 1):
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
                spv[name] =(np.ones((spv['nper'],))*spv[name]).astype(spv[name].dtype)
            else:
                raise ValueError('{} must have the same length as the number of periods ({})'.format(name,spv['nper']))

    return spv


def change_stress_period_settings (m, spv):
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
    spv =_check_stress_period_values(spv)
    dis = m.get_package('dis')
    nper = spv['nper']
    if 'oc_stress_per_data' not in spv.keys():
        oc_stress_per_data = None
    else:
        oc_stress_per_data = spv['oc_stress_per_data']
    dis.nper = nper
    dis.perlen = flopy.utils.Util2d(m,(nper,),np.float32,spv['perlen'],'perlen')
    dis.nstp = flopy.utils.Util2d(m,(nper,),int,spv['nstp'],'nstp')
    dis.steady = flopy.utils.Util2d(m,(nper,),bool,spv['steady'],'steady')
    dis.tsmult = flopy.utils.Util2d(m,(nper,),np.float32,spv['tsmult'],'tsmult')
    dis.check()
    if oc_stress_per_data is None:
        m.oc.stress_period_data = {(0,0):['save head', 'save drawdown', 'save budget', 'print budget']}
    else:
        m.oc.stress_period_data = oc_stress_per_data



def get_base_mf_ss (model_version='m_strong_vert', recalc=False):
    """
    load a steady state modflow model which has been run
    :param recalc:
    :param model_version: which model version to use... see supporting data path for values
    :return:
    """
    base_path = '{}/base_model_runs/{}_base_ss_mf'.format(sdp,model_version)

    if os.path.exists('{}/{}_base_SS.hds'.format(base_path, model_version)) and not recalc:
        m = flopy.modflow.Modflow.load('{}/{}_base_SS.nam'.format(base_path,model_version),model_ws=base_path,forgive=False)
        return m

    m = import_gns_model(name='{}_base_SS'.format(model_version), dir_path=base_path, safe_mode=False, mt3d_link=True,
                         model_version=model_version)
    m.write_input()
    m.write_name_file()
    m.run_model()
    return m


if __name__ == '__main__':
    run_type = 3
    if run_type == 0:
        test_dir = r"C:\Users\MattH\Desktop\test_dir2"
        mf = import_gns_model('test2', test_dir, safe_mode=True)
        mf.write_input()
        mf.write_name_file()
    if run_type == 3:
        get_base_mf_ss()
