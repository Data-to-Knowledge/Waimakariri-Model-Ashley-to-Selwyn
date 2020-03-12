# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/09/2017 8:28 AM
"""

from __future__ import division

import itertools
import os

import numpy as np

import env
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools.metadata_managment.cwms_index import get_zone_array_index
from waimak_extended_boundry.model_run_tools.model_setup.realisation_id import \
    get_rch_multipler
from waimak_extended_boundry.non_model_work.lsr_support.map_rch_to_model_array import \
    map_rch_to_array
from waimak_extended_boundry.model_and_NSMC_build.m_packages.rch_packages import get_rch_fixer, _get_rch

lsrm_rch_base_dir = env.gw_met_data('niwa_netcdf/lsrm/lsrm_results/water_year_means')
rch_idx_shp_path = env.gw_met_data("niwa_netcdf/lsrm/lsrm_results/test/output_test2.shp")


def get_forward_rch(model_id, naturalised, pc5=False, rcm=None, rcp=None, period=None,
                    amag_type=None, cc_to_waimak_only=False, super_gmp=False):
    """
    get the rch for the forward runs
    :param model_id: which NSMC realisation to use
    :param naturalised: boolean if True then get rch for
    :param rcm: regional Climate model identifier
    :param rcp: representetive carbon pathway identifier
    :param period: e.g. 2010, 2020, ect
    :param amag_type: the amalgamation type one of: 'tym': twenty year mean, was 10 but changed on 20/09/2017
                                                    'min': minimum annual average,
                                                    'low_3_m': average of the 3 lowest consecutive years
                                                    'mean': full data mean
                                                    None: then use 'mean'
    :param pc5: boolean if true use assumed PC5 efficency (only applied to the WILS and {something} areas)
    :param cc_to_waimak_only: if true only apply the cc rch to the waimakairi zone (use the vcsn data other wise (keep senario)
    :param super_gmp: boolean if True then use the larger reduction to the will command area
    :return: rch array (11,364,365)
    """
    # I think I need to apply everything as a percent change or somehow normalise the rch so that I do not get any big
    # changes associated with changes in models which created the recharge array.

    hdf_path = _get_rch_hdf_path(lsrm_rch_base_dir, naturalised, pc5, rcm, rcp)
    # get rch array from LSRM

    amalg_dict = {None: 'mean', 'mean': 'mean', 'tym': 'period_mean', 'low_3_m': '3_lowest_con_mean',
                  'min': 'lowest_year'}

    method = amalg_dict[amag_type]  # some of these naming conventions will not owrk now. see _create_all_lsrm_arrays
    sen = 'current'
    if naturalised:
        sen = 'nat'
    if pc5:
        sen = 'pc5'

    rch_array = get_lsrm_base_array(sen, rcp, rcm, period, method)
    if super_gmp:
        assert sen == 'current', 'for super gmp senario must be current'
        assert all([e is None for e in [rcp,
                                        rcm]]), 'rcp, rcm, must all be None, no support for climate change scenarios'
        mult = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\cmp_gmp_point_sources_n.shp".format(smt.sdp),
                                             'drn_change', True)
        mult[np.isnan(mult)]=1

        rch_array *= mult

    # apply multiplier array from pest parameraterisation
    rch_mult = get_rch_multipler(model_id)
    rch_array *= rch_mult

    if cc_to_waimak_only:
        base_rch = _get_rch(2)
        base_rch *= rch_mult
        idx_array = get_zone_array_index(['chch', 'selwyn'])
        rch_array[idx_array] = base_rch[idx_array]
    # handle weirdness from the arrays (e.g. ibound ignore the weirdness from chch/te waihora paw)

    # fix tewai and chch weirdeness
    fixer = get_rch_fixer()
    # chch
    rch_array[fixer == 0] = 0.0002
    # te wai and coastal
    rch_array[fixer == 1] = 0

    no_flow = smt.get_no_flow(0)
    no_flow[no_flow < 0] = 0
    rch_array[~no_flow.astype(bool)] = 0
    return rch_array


def _get_rch_hdf_path(base_dir, naturalised, pc5, rcm, rcp):
    """
    get the path for the rch
    :param base_dir: the directory containing all forward runs
    :param naturalised: boolean if true use array with no irrigation
    :param pc5: boolean if true use 100% effcient irrigation
    :param rcm: None or the RCM of the model
    :param rcp: None or the RCP of the model
    :return:
    """
    if rcp is None and rcm is not None:
        raise ValueError('rcm and rcp must either both be none or be defined')
    elif rcm is None and rcp is not None:
        raise ValueError('rcm and rcp must either both be none or be defined')
    elif rcp is None and rcp is None:
        if naturalised:
            outpath = os.path.join(base_dir, 'wym_vcsn_no_irr.h5')
        elif pc5:
            outpath = os.path.join(base_dir, 'wym_vcsn_100perc.h5')
        else:
            outpath = os.path.join(base_dir, 'wym_vcsn_80perc.h5')
    else:
        if naturalised:
            outpath = os.path.join(base_dir, "wym_{}_{}_no_irr.h5".format(rcp, rcm))
        elif pc5:
            outpath = os.path.join(base_dir, "wym_{}_{}_100perc.h5".format(rcp, rcm))
        else:
            outpath = os.path.join(base_dir, "wym_{}_{}_80perc.h5".format(rcp, rcm))

    return outpath


def _create_all_lsrm_arrays():
    """
    saves arrays for all of the periods and types ect as needed saves both rch and ird arrays
    :return:
    """
    if not os.path.exists(os.path.join(lsrm_rch_base_dir, 'arrays_for_modflow')):
        os.makedirs(os.path.join(lsrm_rch_base_dir, 'arrays_for_modflow'))
    zidx = get_zone_array_index('waimak')
    site_list = list(set(smt.shape_file_to_model_array(rch_idx_shp_path, 'site', True)[zidx]))
    periods = range(2010, 2100, 20)
    rcps = ['RCP4.5', 'RCP8.5']
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    senarios = ['pc5', 'nat', 'current']
    # cc stuff
    for per, rcp, rcm, at, sen in itertools.product(periods, rcps, rcms, amalg_types, senarios):
        print ((per, rcp, rcm, at, sen))
        naturalised = False
        pc5 = False
        if sen == 'nat':
            naturalised = True
        elif sen == 'pc5':
            pc5 = True
        elif sen == 'current':
            pass
        else:
            raise ValueError('shouldnt get here')

        hdf_path = _get_rch_hdf_path(base_dir=lsrm_rch_base_dir, naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp)
        temp, ird = map_rch_to_array(hdf=hdf_path,
                                     method=at,
                                     period_center=per,
                                     mapping_shp=rch_idx_shp_path,
                                     period_length=20,
                                     return_irr_demand=True,
                                     site_list=site_list)
        outpath = os.path.join(lsrm_rch_base_dir,
                               'arrays_for_modflow/rch_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        outpath_ird = os.path.join(lsrm_rch_base_dir,
                                   'arrays_for_modflow/ird_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        np.savetxt(outpath, temp)
        np.savetxt(outpath_ird, ird)

    # RCP past
    amalg_types = ['period_mean', '3_lowest_con_mean', 'lowest_year']
    for rcm, sen, at in itertools.product(rcms, senarios, amalg_types):
        naturalised = False
        pc5 = False
        if sen == 'nat':
            naturalised = True
        elif sen == 'pc5':
            pc5 = True
        elif sen == 'current':
            pass
        else:
            raise ValueError('shouldnt get here')
        per = 1980
        rcp = 'RCPpast'
        print ((per, rcp, rcm, at, sen))
        hdf_path = _get_rch_hdf_path(base_dir=lsrm_rch_base_dir, naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp)
        temp, ird = map_rch_to_array(hdf=hdf_path,
                                     method=at,
                                     period_center=per,
                                     mapping_shp=rch_idx_shp_path,
                                     period_length=50,
                                     return_irr_demand=True)
        outpath = os.path.join(lsrm_rch_base_dir,
                               'arrays_for_modflow/rch_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        outpath_ird = os.path.join(lsrm_rch_base_dir,
                                   'arrays_for_modflow/ird_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        np.savetxt(outpath, temp)
        np.savetxt(outpath_ird, ird)

    # VCSN
    for sen in senarios:
        naturalised = False
        pc5 = False
        if sen == 'nat':
            naturalised = True
        elif sen == 'pc5':
            pc5 = True
        elif sen == 'current':
            pass
        else:
            raise ValueError('shouldnt get here')
        at = 'mean'
        per = None
        rcp = None
        rcm = None
        print ((per, rcp, rcm, at, sen))
        hdf_path = _get_rch_hdf_path(base_dir=lsrm_rch_base_dir, naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp)
        temp, ird = map_rch_to_array(hdf=hdf_path,  # handle IRD as text files for now
                                     method=at,
                                     period_center=per,
                                     mapping_shp=rch_idx_shp_path,
                                     period_length=20,
                                     return_irr_demand=True)
        outpath = os.path.join(lsrm_rch_base_dir,
                               'arrays_for_modflow/rch_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        outpath_ird = os.path.join(lsrm_rch_base_dir,
                                   'arrays_for_modflow/ird_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
        np.savetxt(outpath, temp)
        np.savetxt(outpath_ird, ird)


def get_lsr_base_period_inputs(sen, rcp, rcm, per, at):
    """
    get the LSR comparison period

    :param sen: the senario, senarios = ['pc5', 'nat', 'current']
    :param rcp: rcps = ['RCPpast', 'RCP4.5', 'RCP8.5']
    :param rcm: rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    :param per: None(vcsn), 1980(RCPpast),  periods = range(2010, 2100, 20) (climate change)
    :param at: ['period_mean', '3_lowest_con_mean', 'lowest_year'] (climate change) ['mean'] vcsn
    :return:
    """
    if rcp is None and rcm is None:
        per = None
        at = 'mean'
        sen = 'current'
    elif rcp is not None and rcm is not None:
        rcp = 'RCPpast'
        per = 1980
        at = 'period_mean'  # setting based on the period mean for RCP past for all
    return (sen, rcp, rcm, per, at)


def get_lsrm_base_array(sen, rcp, rcm, per, at):
    """
    get the lsr array
    :param sen: see above
    :param rcp:
    :param rcm:
    :param per:
    :param at:
    :return:
    """
    path = os.path.join(lsrm_rch_base_dir, 'arrays_for_modflow/rch_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
    if not os.path.exists(path):
        raise ValueError('array not implemented, why are you using {}'.format((sen, rcp, rcm, per, at)))
    outdata = np.loadtxt(path)
    if outdata.shape != (smt.rows, smt.cols):
        raise ValueError('incorrect shape for rch array: {}'.format(outdata.shape))

    return outdata


def get_ird_base_array(sen, rcp, rcm, per, at):
    """
    get the irrigation demand array
    :param sen: see above
    :param rcp:
    :param rcm:
    :param per:
    :param at:
    :return:
    """
    path = os.path.join(lsrm_rch_base_dir, 'arrays_for_modflow/ird_{}_{}_{}_{}_{}.txt'.format(sen, rcp, rcm, per, at))
    if not os.path.exists(path):
        raise ValueError('array not implemented, why are you using {}'.format((sen, rcp, rcm, per, at)))
    outdata = np.loadtxt(path)
    if outdata.shape != (smt.rows, smt.cols):
        raise ValueError('incorrect shape for ird: {}'.format(outdata.shape))
    if sen == 'current':
        outdata *= 1.2  # this accounts for the 20 % leakage in our current senario which is 80% efficient.  there is no difference between the two irrigation demand arrays otherwise
    return outdata


if __name__ == '__main__':
    # tests
    testtype = 0
    if testtype == 0:
        waimak = get_zone_array_index('waimak')
        coastal_waimak = get_zone_array_index('coastal_waimak')
        inland_waimak = get_zone_array_index('inland_waimak')
        super_gmp = get_forward_rch(model_id='NsmcBase', naturalised=False, pc5=False, rcm=None, rcp=None, period=None,
                    amag_type=None, cc_to_waimak_only=False, super_gmp=True)
        gmp = get_forward_rch(model_id='NsmcBase', naturalised=False, pc5=True, rcm=None, rcp=None, period=None,
                    amag_type=None, cc_to_waimak_only=False, super_gmp=False)
        cmp_rch = get_forward_rch(model_id='NsmcBase', naturalised=False, pc5=False, rcm=None, rcp=None, period=None,
                    amag_type=None, cc_to_waimak_only=False, super_gmp=False)

        for name in ['waimak', 'inland_waimak', 'coastal_waimak']:
            idx = eval(name)
            print('{} cmp: {} m3/s'.format(name, np.nansum(cmp_rch[idx]*200*200)/86400))
            print('{} super_gmp: {}, m3/s'.format(name, np.nansum(super_gmp[idx]*200*200)/86400))
            print('{} gmp: {}, m3/s'.format(name, np.nansum(gmp[idx]*200*200)/86400))
        print('done')
    if testtype == 1:
        _create_all_lsrm_arrays()

    if testtype == 2:
        {None: 'mean', 'mean': 'mean', 'tym': 'period_mean', 'low_3_m': '3_lowest_con_mean',
         'min': 'lowest_year'}
        periods = range(2010, 2100, 20)
        rcps = ['RCP4.5', 'RCP8.5']
        rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
        amalg_types = ['tym', 'low_3_m', 'min']
        senarios = ['pc5', 'nat', 'current']
        # cc stuff
        for per, rcp, rcm, at, sen in itertools.product(periods, rcps, rcms, amalg_types, senarios):
            naturalised = False
            pc5 = False
            if sen == 'nat':
                naturalised = True
            elif sen == 'pc5':
                pc5 = True
            elif sen == 'current':
                pass

            test = get_forward_rch('opt', naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp, period=per, amag_type=at,
                                   cc_to_waimak_only=True)

        amalg_types = ['tym', 'low_3_m', 'min']
        for rcm, sen, at in itertools.product(rcms, senarios, amalg_types):
            naturalised = False
            pc5 = False
            if sen == 'nat':
                naturalised = True
            elif sen == 'pc5':
                pc5 = True
            elif sen == 'current':
                pass
            else:
                raise ValueError('shouldnt get here')
            per = 1980
            rcp = 'RCPpast'
            get_forward_rch('opt', naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp, period=per, amag_type=at)

        for sen in senarios:
            naturalised = False
            pc5 = False
            if sen == 'nat':
                naturalised = True
            elif sen == 'pc5':
                pc5 = True
            elif sen == 'current':
                pass
            else:
                raise ValueError('shouldnt get here')
            at = 'mean'
            per = None
            rcp = None
            rcm = None
            get_forward_rch('opt', naturalised=naturalised, pc5=pc5, rcm=rcm, rcp=rcp, period=per, amag_type=at)
