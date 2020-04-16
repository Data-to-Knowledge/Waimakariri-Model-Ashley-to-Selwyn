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
import netCDF4 as nc
from waimak_extended_boundry.extended_boundry_model_tools import smt
from waimak_extended_boundry.model_run_tools.metadata_managment.cwms_index import get_zone_array_index
from waimak_extended_boundry.model_run_tools.model_setup.realisation_id import \
    get_rch_multipler

rch_data_path = os.path.join(env.sdp_required,'recharge_arrays.nc')

def get_optimisation_recharge():
    """
    return the base recharge array without any pest multipliers,
    this is the recharge used for the start of the optimisation
    :return:
    """
    data = np.array(nc.Dataset(rch_data_path).variables['opt_rch'])
    return data

def get_rch_fixer():
    """
    an array to index the abnormal recharge in chch and te waihora 1 = tewaihora and coastal, 0 = chch, all others nan
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return:
    """
    data = np.array(nc.Dataset(rch_data_path).variables['recharge_fixer'])
    return data


def get_forward_rch(model_id, naturalised, pc5=False, rcm=None, rcp=None, period=None,
                    amag_type=None, cc_to_waimak_only=False, super_gmp=False):
    """
    get the rch for the forward runs #todo this could actually use more documentation... it's a bit confusing
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
        # this was not extensivly used, i think it was mostly a check
        assert sen == 'current', 'for super gmp senario must be current'
        assert all([e is None for e in [rcp,
                                        rcm]]), 'rcp, rcm, must all be None, no support for climate change scenarios'
        mult = smt.shape_file_to_model_array(os.path.join(env.sdp_required,"shp/cmp_gmp_point_sources_n.shp"),
                                             'drn_change', True)
        mult[np.isnan(mult)]=1

        rch_array *= mult

    # apply multiplier array from pest parameraterisation
    rch_mult = get_rch_multipler(model_id)
    rch_array *= rch_mult

    if cc_to_waimak_only:
        base_rch = get_optimisation_recharge()
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
    get the lsr array see below for requirments
    :param sen:
    :param rcp:
    :param rcm:
    :param per:
    :param at:
    :return:
    """
    return _get_rch_ird(sen,rcp,rcm,per,at,recharge=True)



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
    return _get_rch_ird(sen,rcp,rcm,per,at,recharge=False)


def _get_rch_ird(sen, rcp, rcm, per, at, recharge):
    """

    :param sen:
    :param rcp:
    :param rcm:
    :param per:
    :param at:
    :param recharge:boolean if True then Recharge else IRD
    :return:
    """
    if recharge:
        key_val = 'recharge'
    else:
        key_val = 'ird'

    senarios = ['pc5', 'nat', 'current']
    rcps = ['RCPpast', 'RCP4.5', 'RCP8.5', None]
    rcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M', None]
    periods = [None, 1980] + range(2010, 2100, 20)
    ats = ['period_mean', '3_lowest_con_mean', 'lowest_year', 'mean']

    # check arguments
    assert sen in senarios
    assert rcp in rcps
    assert rcm in rcms
    assert per in periods
    assert at in ats

    data = nc.Dataset(rch_data_path)

    # capture the current data
    if rcp is None and rcm is None:
        assert per is None
        assert at is 'mean'
        idx = np.where(np.array(data['scenario']) == sen)[0][0]
        out = np.array(np.array(data['current_{}'.format(key_val)])[idx])

    # capture the rcp_past
    elif rcp == 'RCPpast':
        assert per == 1980
        assert at != 'mean'
        i = np.where(np.array(data.variables['scenario']) == sen)[0][0]
        j = np.where(np.array(data.variables['rcm']) == rcm)[0][0]
        k = np.where(np.array(data.variables['amalg_type']) == at)[0][0]
        out = np.array(data.variables['rcp_past_{}'.format(key_val)][i, j, k])

    # climate change scenarios
    elif rcp is not None and rcm is not None:
        assert per in range(2010, 2100, 20)
        assert at != 'mean'

        i = np.where(np.array(data.variables['scenario']) == sen)[0][0]
        j = np.where(np.array(data.variables['rcp']) == rcp)[0][0]
        k = np.where(np.array(data.variables['rcm']) == rcm)[0][0]
        l = np.where(np.array(data.variables['period']) == per)[0][0]
        m = np.where(np.array(data.variables['amalg_type']) == at)[0][0]
        out = np.array(data.variables['future_{}'.format(key_val)][i, j, k, l, m])

    else:
        raise ValueError('rcm and rcp must either be both None or both not None')

    assert out.shape == (smt.rows, smt.cols), 'weird shape {}'.format(out.shape)

    return out

# deprecidated functions and opperations

rch_idx_shp_path = env.gw_met_data("niwa_netcdf/lsrm/lsrm_results/test/output_test2.shp")

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
    raise NotImplementedError('left for documentation purposes only')
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
    raise NotImplementedError('this is left of documentation purposes only')
    from waimak_extended_boundry.non_model_work.lsr_support.map_rch_to_model_array import \
        map_rch_to_array
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


if __name__ == '__main__':
    # tests
    pass

