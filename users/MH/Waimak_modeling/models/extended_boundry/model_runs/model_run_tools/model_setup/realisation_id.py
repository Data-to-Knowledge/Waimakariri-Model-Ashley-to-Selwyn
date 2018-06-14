# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/09/2017 2:23 PM
"""

from __future__ import division
from core import env
import flopy_mh as flopy
import pickle
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import os
from warnings import warn
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import _get_wel_spd_v1, _get_wel_spd_v2
import pandas as pd
import numpy as np
import subprocess
import netCDF4 as nc
import shutil
from users.MH.Waimak_modeling.supporting_data_path import sdp
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.convergance_check import \
    modflow_converged
from copy import deepcopy

temp_pickle_dir = '{}/temp_pickle_dir'.format(smt.pickle_dir)
if not os.path.exists(temp_pickle_dir):
    os.makedirs(temp_pickle_dir)


def get_base_well(model_id, org_pumping_wells, recalc=False):
    """
    applies the NSMC pumping mulitpliers
    :param model_id: the NSMC realisation
    :param org_pumping_wells: if True use the model peiod wells if false use the 2014-2015 usage
    :param recalc: usual recalc
    :return:
    """
    if org_pumping_wells:
        org = 'model_period_pumping'
    else:
        org = '14-15_pumping_waimak'
    pickle_path = "{}/model_{}_base_wells_{}.p".format(temp_pickle_dir, model_id, org)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata

    if org_pumping_wells:
        all_wells = _get_wel_spd_v1()  # usage for the model period as passed to the optimisation
    else:
        all_wells = _get_wel_spd_v2()  # usage for 2014/2015 period in waimak zone

    all_wells.loc[:, 'nsmc_type'] = ''

    # pumping wells
    all_wells.loc[(all_wells.type == 'well') & (all_wells.cwms == 'chch'), 'nsmc_type'] = 'pump_c'
    all_wells.loc[(all_wells.type == 'well') & (all_wells.cwms == 'selwyn'), 'nsmc_type'] = 'pump_s'
    all_wells.loc[(all_wells.type == 'well') & (all_wells.cwms == 'waimak'), 'nsmc_type'] = 'pump_w'

    # selwyn_hillfeds
    all_wells.loc[all_wells.type == 'river', 'nsmc_type'] = 'sriv'

    # races
    all_wells.loc[(all_wells.type == 'race') & (all_wells.zone == 'n_wai'), 'nsmc_type'] = 'n_race'
    all_wells.loc[(all_wells.type == 'race') & (all_wells.zone == 's_wai'), 'nsmc_type'] = 's_race'

    # boundary fluxes
    all_wells.loc[all_wells.type == 'llr_boundry_flux', 'nsmc_type'] = 'llrzf'
    all_wells.loc[all_wells.type == 'ulr_boundry_flux', 'nsmc_type'] = 'ulrzf'
    all_wells.loc[all_wells.type == 'boundry_flux', 'nsmc_type'] = 'nbndf'

    base_dir = os.path.dirname(get_model_name_path(model_id))
    mult_path = '{}/wel_adj.txt'.format(base_dir)  # figure out what the pest control file will populate
    multipliers = pd.read_table(mult_path, index_col=0, delim_whitespace=True, names=['value'])

    mult_groups = ['pump_c', 'pump_s', 'pump_w', 'sriv', 'n_race', 's_race', 'nbndf']
    for group in mult_groups:
        all_wells.loc[all_wells.nsmc_type == group, 'flux'] *= multipliers.loc[group, 'value']

    add_groups = ['llrzf', 'ulrzf']
    for group in add_groups:
        all_wells.loc[all_wells.nsmc_type == group, 'flux'] = multipliers.loc[group, 'value'] / (
            all_wells.nsmc_type == group).sum()
    pickle.dump(all_wells, open(pickle_path, 'w'))
    return all_wells


def get_base_rch(model_id, recalc=False):
    """
    loads the rch from the model and saves as a pickle for faster use # depreciated
    :param model_id: the model realisation to use
    :param recalc: False
    :return:
    """
    pickle_path = "{}/model_{}_base_rch.p".format(temp_pickle_dir, model_id)
    if (os.path.exists(pickle_path)) and (not recalc):
        rch = pickle.load(open(pickle_path))
        return rch

    m = get_model(model_id)
    rch = m.rch.rech.array[0, 0, :, :]
    pickle.dump(rch, open(pickle_path, 'w'))
    return rch


def get_rch_multipler(model_id):
    """
    get the recharge multipler if it does not exist in the file then create it with fac2real
    :param model_id: the NSMC realisation
    :return:
    """
    model_ws = os.path.dirname(get_model_name_path(model_id))
    rch_mult_path = os.path.join(model_ws, 'recharge_mul.ref')
    if not os.path.exists(rch_mult_path):
        model_ws = os.path.splitunc(model_ws)[1].replace('/SCI',
                                                         'P:/')  # for now assuming that SCI is mapped to P drive could fix in future
        exe_path = r"P:/Groundwater/Matt_Hanson/model_exes/fac2real.exe"
        test = subprocess.call(exe_path + ' < fac2real_rech.in', cwd=model_ws, shell=True)
        print(test)

    # to read in the possibly wrapped array use:
    outdata = flopy.utils.Util2d.load_txt((smt.rows, smt.cols), rch_mult_path, float, '(FREE)')
    no_flow = smt.get_no_flow(1)
    no_flow[no_flow < 0] = 0
    outdata[~no_flow.astype(bool)] = 1
    return outdata


def get_model_name_path(model_id):
    """
    get the path to a model_id base model
    # !!!model id needs to be non-numeric (start with a letter) and does not contain an '_' !!!!
    :param model_id:
    :return:
    """
    # new model check the well package, check the sfr package and run new_model_run

    if model_id in ['test', 'opt', 'NsmcBaseB', 'VertUnstabA', 'VertUnstabB']:
        warn('model {} is depreciated'.format(model_id))

    model_dict = {
        # a place holder to test the scripts
        'test': r"C:\Users\MattH\Desktop\Waimak_modeling\ex_bd_tester\test_import_gns_mod\mf_aw_ex.nam",

        # the optimized model as of 26/09/2017; depreciated due to concerns in eyrewell forrest area and changes
        # in the recharge and well packages in the waimakariri Zone
        'opt': "{}/from_gns/optimised/mf_aw_ex/mf_aw_ex.nam".format(smt.sdp),

        # an opitimised model recieved 20/10/2017 that has the new priors, well and rch packages
        # which fits the heads, streams, ofshore flows, but not the deep heads, and not the verticle targets
        'StrOpt': "{}/from_gns/StrOpt/AW20171019_i3_optver/mf_aw_ex.nam".format(smt.sdp),

        # an optimised model revieved 24/10/2017 which fits all targets, but is rather unstable
        'VertUnstabB': "{}/from_gns/VertUnstabB/AW20171022_i4_optver/i4/mf_aw_ex.nam".format(smt.sdp),

        # the previous optimisation recived 24/10/2017 iterattion to VertUnstabB which does not fit the targets as well,
        # but is more stable
        'VertUnstabA': "{}/from_gns/VertUnstabA/AW20171022_i3_optver/i3/mf_aw_ex.nam".format(smt.sdp),

        # optimisation recieved on 24/10/2017 which is more stable and better hits the targets than either of the
        # VertUnstab models.  if the jacobian runs fine then this will be the base for the NSMC this ended up the base for the NSMC
        'NsmcBase': "{}/from_gns/NsmcBase/AW20171024_2_i2_optver/i2/mf_aw_ex.nam".format(smt.sdp),

        # optimisation revieved on 25/10/2017 which is the previous iteration to NsmcBase.  NsmcBase was quite unstable,
        # which is why we are considering the previous optimisation
        'NsmcBaseB': "{}/from_gns/NsmcBaseB/AW20171024_2_i1_optver/i1/mf_aw_ex.nam".format(smt.sdp),

        # optimisation recieved on 8/01/2018 which was run over xmas to address concerns over the amount of flow lost
        # in the upper ashley basin e.g. near oxford the model may be overfit and it is quite unstable as such is has
        # been suggested as a sensitivity analysis for this concern
        'AshOpt': "{}/from_gns/AshOpt/AW20180103_Ash0_Opt/AW20180103_Ash0_Opt/mf_aw_ex.nam".format(smt.sdp)
    }
    if '_' in model_id:
        raise ValueError(
            '_ in model id: {}, model_id cannot include an "_" as this is a splitting character'.format(model_id))
    if model_id not in model_dict.keys() and not 'NsmcReal' in model_id:
        raise NotImplementedError('model {} has not yet been defined'.format(model_id))

    if 'NsmcReal' in model_id:
        m = get_model(model_id, save_to_dir=True)
        name_path = m.namefile
    else:
        name_path = model_dict[model_id]
    return name_path


def _get_nsmc_realisation(model_id, save_to_dir=False):
    """
    wrapper to get model from a NSMC realisation
    :param model_id: identifier 'NsmcReal{nsmc_num:06d}'
    :param save_to_dir: boolean if true save a copy of the model for quicker reteval in the dir specified below
    :param temp_int: a temp interger to make unique temporaty files only needed to handle multiprocessing applications
    :return:
    """
    assert 'NsmcReal' in model_id, 'unknown model id: {}, expected NsmcReal(nsmc_num:06d)'.format(model_id)
    assert len(model_id) == 14, 'unknown model id: {}, expected NsmcReal(nsmc_num:06d)'.format(model_id)
    base_converter_dir = "{}/base_for_nsmc_real".format(smt.sdp)
    # check if the model has previously been saved to the save dir, and if so, load from there
    save_dir = env.gw_met_data("mh_modeling/nsmc_loaded_realisations_TEMP")
    converter_dir = os.path.join(os.path.expanduser('~'), 'temp_nsmc_generation{}'.format(os.getpid()))
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if os.path.exists(os.path.join(save_dir, '{}_base'.format(model_id),'{}_base.hds'.format(model_id))):
        name_file_path = os.path.join(save_dir, '{}_base'.format(model_id), '{m}_base.nam'.format(m=model_id))
        m = flopy.modflow.Modflow.load(name_file_path, model_ws=os.path.dirname(name_file_path), forgive=False)
        return m

    # copy the orginal converter dir to the temporary working dir
    shutil.copytree(base_converter_dir, converter_dir)

    nsmc_num = int(model_id[-6:])
    param_data = nc.Dataset(env.gw_met_data(r"mh_modeling\netcdfs_of_key_modeling_data\nsmc_params_obs_metadata.nc"))
    param_idx = np.where(np.array(param_data.variables['nsmc_num']) == nsmc_num)[0][0]

    print('writing data to parameter files')
    # write well paraemeters to wel_adj.txt
    with open(os.path.join(converter_dir, 'wel_adj.txt'), 'w') as f:
        for param in ['pump_c', 'pump_s', 'pump_w', 'sriv', 'n_race', 's_race', 'nbndf', 'llrzf', 'ulrzf']:
            val = param_data.variables[param][param_idx]
            f.write('{} {}\n'.format(param, val))

    # write rch parameters to rch_ppts.txt
    with open(os.path.join(converter_dir, 'rch_ppts.txt'), 'w') as f:
        keys = np.array(param_data.variables['rch_ppt'])
        x = np.array(param_data.variables['rch_ppt_x'])
        y = np.array(param_data.variables['rch_ppt_y'])
        group = np.array(param_data.variables['rch_ppt_group'])
        val = np.array(param_data.variables['rch_mult'][param_idx])
        for k, _x, _y, g, v in zip(keys, x, y, group, val):
            f.write('{} {} {} {} {}\n'.format(k, _x, _y, g, v))

    # write kh and kv
    all_kv = np.array(param_data.variables['kv'][param_idx])  # shape = (11, 178)
    all_kh = np.array(param_data.variables['kh'][param_idx])  # shape = (11, 178)
    all_names = np.array(param_data.variables['khv_ppt'])
    all_x = np.array(param_data.variables['khv_ppt_x'])
    all_y = np.array(param_data.variables['khv_ppt_y'])
    for layer in range(smt.layers):
        khv_idx = np.isfinite(all_kh[layer])
        layer_names = all_names[khv_idx]
        layer_x = all_x[khv_idx]
        layer_y = all_y[khv_idx]
        layer_kv = all_kv[layer][khv_idx]
        layer_kh = all_kh[layer][khv_idx]

        # parameters to kh_kkp_{layer one indexed}.txt
        with open(os.path.join(converter_dir, 'KH_ppk_{:02d}.txt'.format(layer + 1)), 'w') as f:  # one indexed
            for n, x, y, kh in zip(layer_names, layer_x, layer_y, layer_kh):
                f.write('{}\t{}\t{}\t1\t{}\n'.format(n, x, y, kh))

        # write kv parameters to kv_ppk_{layer one indexed}.txt
        with open(os.path.join(converter_dir, 'KV_ppk_{:02d}.txt'.format(layer + 1)), 'w') as f:  # one indexed
            for n, x, y, kv in zip(layer_names, layer_x, layer_y, layer_kv):
                f.write('{}\t{}\t{}\t1\t{}\n'.format(n, x, y, kv))

    # write sfr parameters to segfile.txt and/or segfile.tpl
    # load template
    sfr_template = pd.read_table(os.path.join(converter_dir, 'sfr_segdata.tpl'), skiprows=1, sep='\t')
    # make dictionary of parameters
    flows = ['$mid_c_flo $', '$top_c_flo $', '$top_e_flo $']
    replacement = {}
    # flows
    for f in flows:
        replacement[f] = param_data.variables[f.strip('$ ')][param_idx]
    # ks
    names = param_data.variables['sfr_cond'][:]
    for k in set(sfr_template.hcond2).union(set(sfr_template.hcond1)):
        name_idx = np.where(names == k.strip(' $').lower())[0][0]
        val = param_data.variables['sfr_cond_val'][param_idx][name_idx]
        replacement[k] = val

    # do the replacement
    sfr_out = sfr_template.replace(replacement)
    sfr_out.to_csv(os.path.join(converter_dir, 'sfr_segdata.txt'), sep='\t', index=False)

    # write fault parameters to fault_ks.txt
    with open(os.path.join(converter_dir, 'fault_ks.txt'), 'w') as f:
        for param in ['fkh_mult', 'fkv_mult']:
            val = param_data.variables[param][param_idx]
            f.write('{}\n'.format(val))

    # write drain package from parameters and mf_aw_ex_drn.tpl

    # make replacement dictionary
    names = {'$   d_ash_c$', '$   d_ash_s$', '$  d_chch_c$', '$  d_cust_c$', '$  d_dlin_c$', '$  d_dsel_c$',
             '$  d_smiths$', '$  d_ulin_c$', '$  d_usel_c$', '$ d_ash_est$', '$ d_cam_yng$', '$ d_dwaimak$',
             '$ d_kairaki$', '$ d_tar_gre$', '$ d_uwaimak$', '$ d_waihora$', '$d_bul_avon$', '$d_bul_styx$',
             '$d_cam_mrsh$', '$d_cam_revl$', '$d_cour_nrd$', '$d_emd_gard$', '$d_kuku_leg$', '$d_nbk_mrsh$',
             '$d_oho_btch$', '$d_oho_jefs$', '$d_oho_kpoi$', '$d_oho_misc$', '$d_oho_mlbk$', '$d_oho_whit$',
             '$d_salt_fct$', '$d_salt_top$', '$d_sbk_mrsh$', '$d_sil_harp$', '$d_sil_heyw$', '$d_sil_ilnd$',
             '$d_tar_stok$'}

    param_names = param_data.variables['drns'][:]
    replacer = {}
    for nm in names:
        drn_idx = param_names == nm.strip('$ ')
        replacer[nm] = param_data.variables['drn_cond'][param_idx][drn_idx][0]

    # replace data
    drn_tpl_path = os.path.join(converter_dir, "mf_aw_ex_drn.tpl")
    with open(drn_tpl_path) as f:
        drns = f.read()
        for key, val in replacer.items():
            drns = drns.replace(key,str(val))
    drns = drns.replace('ptf $\n','') # get rid of the header for pest

    # write new data
    out_drn_path = os.path.join(converter_dir, 'mf_aw_ex.drn')
    with open(out_drn_path,'w') as f:
        f.write(drns)


    # run model.bat
    print('running model.bat')
    cwd = os.path.splitunc(converter_dir)[1].replace('/SCI', 'P:/')
    # for now assuming that SCI is mapped to P drive could fix in future

    p = subprocess.Popen([os.path.join(cwd, "model.bat")], cwd=cwd, shell=True)
    out = p.communicate()
    print(out)
    p = subprocess.Popen(['python', os.path.join(cwd, "aw_gen_faultreal.py")], cwd=cwd)
    p.communicate()
    p = subprocess.Popen(['python', os.path.join(cwd, "aw_gen_sfr.py")], cwd=cwd)
    p.communicate()
    p = subprocess.Popen(['python', os.path.join(cwd, "well_pkg_adjust.py")], cwd=cwd)
    p.communicate()
    # load model
    name_file_path = os.path.join(converter_dir, 'mf_aw_ex.nam')
    # this is a bit unnecessary, but I want it to be exactly the same as the other loader
    m = flopy.modflow.Modflow.load(name_file_path, model_ws=os.path.dirname(name_file_path), forgive=False)

    # if save then save to a new folder and run the model
    if save_to_dir:
        name = '{}_base'.format(model_id)
        dir_path = os.path.join(save_dir, name)
        print('saving model to holding dir'.format(dir_path))
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)  # remove old files to prevent file mix ups
        m._set_name(name)
        m.exe_name = "{}/models_exes/MODFLOW-NWT_1.1.2/MODFLOW-NWT_1.1.2/bin/MODFLOW-NWT_64.exe".format(sdp)
        units = deepcopy(m.output_units)
        for u in units:
            m.remove_output(unit=u)

        fnames = [m.name + e for e in ['.hds', '.ddn', '.cbc', '.sfo']]  # output extension
        funits = [30, 31, 740, 741]  # fortran unit
        fbinflag = [True, True, True, True, ]  # is binary
        fpackage = [[], [], ['UPW', 'DRN', 'RCH', 'SFR', 'WEL'], ['SFR']]
        for fn, fu, fb, fp, in zip(fnames, funits, fbinflag, fpackage):
            m.add_output(fn, fu, fb, fp)

        m.write_name_file()
        m.write_input()
        success, buff = m.run_model()
        con = modflow_converged(os.path.join(converter_dir, m.lst.file_name[0]))
        m.change_model_ws(dir_path) #todo flopy update add an update that changes teh namefile when the directory changes
        m.namefile = os.path.join(m.model_ws, m.namefile)
        if not con or not success:
            os.remove(os.path.join(converter_dir, '{}.hds'.format(m.name)))
            shutil.rmtree(converter_dir)
            raise ValueError('the model did not converge: \n'
                             '{}\n, headfile deleted to prevent running'.format(os.path.join(dir_path, name)))
        shutil.copytree(converter_dir, dir_path)
    shutil.rmtree(converter_dir)
    return m


def get_model(model_id, save_to_dir=False):
    """
    load a flopy model instance of the model id
    :param model_id:
    :param save_to_dir: boolean only for nsmc realisations save the model to dir?
    :return: m flopy model instance
    """
    # check well packaged loads appropriately! yep this is a problem because we define the iface value,
    #  but it isn't used, remove iface manually
    # if it doesn't exist I will need to add the options flags to SFR package manually

    if 'NsmcReal' in model_id:  # model id should be 06d
        assert len(model_id) == 14, 'model id for realsiation is "NsmcReal{nsmc_num:06d}"'
        m = _get_nsmc_realisation(model_id, save_to_dir)
        return m

    name_file_path = get_model_name_path(model_id)
    well_path = name_file_path.replace('.nam', '.wel')
    with open(well_path) as f:
        counter = 0
        while counter < 10:
            line = f.readline()
            counter += 1
            if 'aux' in line.lower():
                raise ValueError(
                    'AUX in well package will cause reading error {} please remove manually'.format(well_path))

    # check SFR package is correct
    sfr_path = name_file_path.replace('.nam', '.sfr')
    with open(sfr_path) as f:
        lines = [next(f).lower().strip() for x in range(10)]

    if any(~np.in1d(['options', 'reachinput', 'end'], lines)):
        raise ValueError('options needed in sfr package {}, please add manually'.format(sfr_path))

    m = flopy.modflow.Modflow.load(name_file_path, model_ws=os.path.dirname(name_file_path), forgive=False)
    return m

def get_stocastic_set(return_model_ids=True):
    """
    a quick wrapper to get all of the 165 stocastic models
    :return:
    """
    nsmc_nums = [5, 17, 18, 26, 37, 44, 72, 103, 117, 133, 142,
                 204, 233, 240, 258, 271, 278, 308, 314, 388, 391, 439,
                 441, 447, 488, 491, 552, 555, 604, 640, 666, 790, 794,
                 800, 808, 809, 871, 883, 932, 952, 955, 992, 1001, 1029,
                 1052, 1055, 1078, 1099, 1106, 1111, 1112, 1114, 1129, 1168, 1181,
                 1250, 1326, 1328, 1330, 1349, 1365, 1368, 1370, 1395, 1399, 1454,
                 1481, 1489, 1498, 1499, 1507, 1574, 1585, 1595, 1618, 1636, 1656,
                 1660, 1681, 1694, 1749, 1773, 1817, 1819, 1836, 1842, 1851, 1870,
                 1881, 1897, 1915, 1916, 1967, 2044, 2052, 2068, 2085, 2114, 2117,
                 2139, 2160, 2200, 2209, 2234, 2278, 2382, 2390, 2430, 2460, 2520,
                 2579, 2588, 2590, 2604, 2643, 2654, 2661, 2673, 2691, 2735, 2748,
                 2752, 2757, 2760, 2818, 2840, 2871, 2875, 2885, 2918, 2971, 2980,
                 2997, 3039, 3064, 3067, 3075, 3100, 3116, 3165, 3216, 3289, 3310,
                 3318, 3350, 3378, 3389, 3418, 3428, 3456, 3464, 3482, 3513, 3585,
                 3641, 3651, 3663, 3685, 3712, 3717, 3762, 3796, 3836, 3861, 3910]

    model_ids = ['NsmcReal{:06d}'.format(e) for e in nsmc_nums]
    if return_model_ids:
        return model_ids
    else:
        return nsmc_nums
if __name__ == '__main__':
    # tests
    m=get_model('NsmcBase')
    print(m)
    print('done')
