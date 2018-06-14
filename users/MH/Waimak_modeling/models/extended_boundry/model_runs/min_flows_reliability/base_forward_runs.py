# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/09/2017 8:39 AM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools import mod_gns_model, \
    zip_non_essential_files, import_gns_model
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools import get_forward_wells, \
    get_forward_rch
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import flopy_mh as flopy
import os
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.convergance_check import modflow_converged
from traceback import format_exc
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.stream_depletion_assesment.raising_heads_no_carpet import \
    get_drn_no_ncarpet_spd
import pandas as pd


def setup_run_forward_run_mp(kwargs):
    """
     a wrapper to allow easy multiprocessing
    :param kwargs: kwargs to pass to setup_run_forward_model
    :return: (model_name, (convergence('convereged'/'did not converge') or trace back to exception))
    """
    try:
        name, success = setup_run_forward_run(**kwargs)
    except Exception as val:
        name = '{}_{}'.format(kwargs['model_id'], kwargs['name'])
        success = format_exc().replace('\n', '')
    return name, success


def setup_run_forward_run(model_id, name, base_dir, cc_inputs=None, pc5=False, pc5_well_reduction=False,
                          pc5_to_waimak_only=False, wil_eff=1, naturalised=False,
                          full_abs=False, pumping_well_scale=1, full_allo=False, org_efficency=None,
                          org_pumping_wells=False, rm_ncarpet=True, super_gmp=False, write_ftl=False, increase_eyre=0):
    """
    sets up and runs a forward run with a number of options
    :param model_id: which NSMC version to user (see mod_gns_model)
    :param name: the name for the model (see mod_gns_model)
    :param base_dir: the folder for the model (see mod_gns_modle)
    :param cc_inputs: a dictionary of:  rcm: regional Climate model identifier,
                                        rcp: representetive carbon pathway identifier,
                                        period: e.g. 2010, 2020, ect,
                                        amag_type: the amalgamation type one of:
                                                        tym: twenty year mean,
                                                        min: minimum annual average,
                                                        low_3_m: average of the 3 lowest consecutive years
                                        cc_to_waimak_only: if True only apply the cc changes to waimak zone
                                                           (still apply senarios e.g. pc5)
    :param pc5: boolean if true use assumed PC5 efficency
    :param pc5_well_reduction: boolean if True, pumping wells in waimakariri will be reduced by 25%
    :param pc5_to_waimak_only: boolean if True, pc5 is only applied to the waimakariri zone
    :param wil_eff: a factor (>0) which will be applied to the WIL scheme losses (e.g. .2 = 20% of optimised model value)
    :param naturalised: boolean if true run with a naturalised system no pumping, irrigation LSR impacts, WIL losses etc.
    :param full_abs: boolean, if true include full abstration for the wells
    :param pumping_well_scale: as factor that is applied to the pumping wells (e.g. .2 = 20%)
                               this factor is separate to any scaling due to additional demand under climate
                               change senarios (this is handled when getting the well data)
                               this factor applies to all wells
    :param full_allo: boolean if true use the full allocation of pumping
    :param org_efficency: not used, held to prevent cleaning up!
    :param org_pumping_wells: if True use the model peiod wells if false use the 2014-2015 usage for the waimak wells
    :param rm_ncarpet: boolean if True remove the N carpet drains
    :param super_gmp: boolean if True run further GMP reductions for will area
    :param write_ftl: boolean if True the model will write the FTL
    :param increase_eyre: a value to increase the eyre inflow in cumics
    :return: (model name, convergence('convereged'/'did not converge'))
    """

    # cc inputs are a dict
    if cc_inputs is None:
        cc_inputs = {'rcm': None, 'rcp': None, 'period': None, 'amag_type': None}
    if not isinstance(cc_inputs, dict):
        raise ValueError('incorrect type for cc_inputs {} expected dict or None'.format(type(cc_inputs)))

    if pc5_to_waimak_only and not pc5:
        raise ValueError('pc5 must be True if pc5 to waimak only is True')
    if pc5_well_reduction and not pc5:
        raise ValueError('pc5 must be True if pc5_well_reduction is True')

    if pc5_to_waimak_only:
        cc_inputs['cc_to_waimak_only'] = True
    well_data, cc_mult, new_water = get_forward_wells(model_id=model_id, full_abstraction=full_abs,
                                                      cc_inputs=cc_inputs, naturalised=naturalised, full_allo=full_allo,
                                                      pc5=pc5_well_reduction, org_pumping_wells=org_pumping_wells)
    well_data.loc[well_data.type == 'well', 'flux'] *= pumping_well_scale
    well_data.loc[(well_data.type == 'race') & (well_data.zone == 'n_wai'), 'flux'] *= wil_eff
    well_data = smt.convert_well_data_to_stresspd(well_data)

    rch = get_forward_rch(model_id, naturalised, pc5, super_gmp=super_gmp, **cc_inputs)

    # I'm assuming that the stream package will not change
    if rm_ncarpet:
        drn = {0: get_drn_no_ncarpet_spd(model_id)}
    else:
        drn = None
    if 'mod_period' in name:
        m = mod_gns_model(model_id,
                          name,
                          base_dir,
                          drain=drn,
                          safe_mode=False,
                          mt3d_link=write_ftl,
                          set_hdry=False)
    else:
        m = mod_gns_model(model_id,
                          name,
                          base_dir,
                          well={0: well_data},
                          recharge={0: rch},
                          drain=drn,
                          safe_mode=False,
                          mt3d_link=write_ftl,
                          set_hdry=False)

    # below included for easy manipulation if needed
    flopy.modflow.mfnwt.ModflowNwt(m,
                                   headtol=1e-5,
                                   fluxtol=500,
                                   maxiterout=100,
                                   thickfact=1e-05,
                                   linmeth=1,
                                   iprnwt=0,  # could change
                                   ibotav=0,
                                   options='COMPLEX',
                                   Continue=False,  # could change
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

    with open(os.path.join(m.model_ws, 'cc_mult_info.txt'), 'w') as f:
        f.write('cc_mult (unitless): ' + str(cc_mult) + '\n')
        f.write('missing_water (m3/day): ' + str(new_water) + '\n')

    # add the extra Eyre water at start of segment 4 near oxford township
    seg_data = pd.DataFrame(m.sfr.segment_data[0])
    seg_data.loc[seg_data.nseg==4,'flow'] += increase_eyre*86400
    m.sfr.segment_data[0] = seg_data.to_records(index=False)

    m.write_name_file()
    m.write_input()
    success, buff = m.run_model(report=True)
    con = False
    if success:
        con = modflow_converged(os.path.join(m.model_ws, m.namefile.replace('.nam', '.list')))
        zip_non_essential_files(m.model_ws, include_list=True)
    if con is None:
        success = 'convergence unknown'
    elif con:
        success = 'converged'
    else:
        success = 'did not converge'
    return '{}_{}'.format(model_id,name), success
