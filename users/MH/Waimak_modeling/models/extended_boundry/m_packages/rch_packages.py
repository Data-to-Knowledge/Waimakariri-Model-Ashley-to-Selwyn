# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/06/2017 11:59 AM
"""

from __future__ import division

import os
import pickle
from copy import deepcopy

import flopy
import numpy as np

from users.MH.Waimak_modeling.model_tools import get_base_rch, no_flow as old_no_flow
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.lsr_support.map_rch_to_model_array import \
    map_rch_to_array
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.lsr_support.generate_an_mean_rch import \
    gen_water_year_average_lsr_irr
import pandas as pd
from warnings import warn


def create_rch_package(m):
    """
    create and add the recharge package
    :param m: a flopy model instance
    :return:
    """
    rch = flopy.modflow.mfrch.ModflowRch(m,
                                         nrchop=3,
                                         ipakcb=740,
                                         rech=_get_rch(),
                                         irch=0,
                                         unitnumber=716)


def _get_rch(version=1, recalc=False):
    """
    a wrapper to get the correct recharge
    :param version: which version to use
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return:
    """
    if version == 1:
        out_rch = _get_rch_v1(recalc)
    elif version == 2:
        out_rch = _get_rch_v2(recalc)
    else:
        raise NotImplementedError('version {} not implemented'.format(version))
    return out_rch


def _get_rch_v1(recalc=False):
    """
    a meld of fouad's LSR mike shear model and expert judgments
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return: rch array m/day (i,j)
    """
    warn('v1 rch is depreciated in newest version of optimised model')
    pickle_path = '{}/org_rch.p'.format(smt.pickle_dir)
    if os.path.exists(pickle_path) and not recalc:
        rch = pickle.load(open(pickle_path))
        return rch

    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = 0
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    confin_rch_zone = smt.shape_file_to_model_array(
        "{}/m_ex_bd_inputs/shp/chch_wm_rch_split_chch_form.shp".format(smt.sdp), 'ID', True)
    zones[(zones == 7) & (np.isfinite(confin_rch_zone))] = 9
    part_zones = deepcopy(zones[:190, :])
    part_zones[old_no_flow] = np.nan  # here this is the no flow for layer 1 I created hence why i don't invert it

    old_rch = get_base_rch()
    scaled_old_rch = np.zeros((190, 365)) * np.nan
    scaled_old_rch[part_zones == 4] = old_rch[part_zones == 4]  # do not scale waimak
    new_old_top = np.nanpercentile(scaled_old_rch, 99)
    scaled_old_rch[scaled_old_rch > new_old_top] = new_old_top
    # create rch values for south wai
    """use homogeneous rate of 270 mm/year based on Williams 2010 (modified from White 2008)
    estimate of 23.8 m³/s for 276,000 ha Te Waihora catchment area.
    Use 190 mm/year for Christchurch WM zone. scaled by an implementation of the david scott model"""

    ds_rch = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp), 'rch',
                                           True) / 1000
    ds_paw = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/dave_scott_rch.shp".format(smt.sdp), 'paw', True)
    # do a bit of cleaning for the DS model
    ds_rch[ds_paw <= 600000] = np.nan  # get rid or rch in streams, te wai, and chch urban area
    new_top = np.nanpercentile(ds_rch, 99)
    ds_rch[ds_rch > new_top] = new_top
    new_bot = np.nanpercentile(ds_rch, 1)
    ds_rch[ds_rch < new_bot] = new_bot
    ds_rch[np.isnan(ds_rch)] = new_bot

    rch = np.zeros((smt.rows, smt.cols))
    rch[zones == 7] = ds_rch[zones == 7] * 175 / 1000 / 365 / ds_rch[zones == 7].mean()
    rch[zones == 9] = ds_rch[zones == 9] * 100 / 1000 / 365 / ds_rch[zones == 9].mean()
    rch[zones == 8] = ds_rch[zones == 8] * 195 / 1000 / 365 / ds_rch[zones == 8].mean()
    rch[zones == 4] = ds_rch[zones == 4] * 290 / 1000 / 365 / ds_rch[zones == 4].mean()
    idx = np.where(np.isfinite(scaled_old_rch))
    rch[idx] = scaled_old_rch[idx]
    # get new rch values for nwai
    pickle.dump(rch, open(pickle_path, 'w'))
    return rch


def _get_rch_v2(recalc=False):
    """
    dave scott recharge for all using 2012 irrigation layer (with some simplifications
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return: rch array m/day (i,j)
    """
    pickle_path = '{}/org_rch_v2.p'.format(smt.pickle_dir)
    if os.path.exists(pickle_path) and not recalc:
        rch = pickle.load(open(pickle_path))
        return rch

    new_no_flow = smt.get_no_flow()
    path = "{}/m_ex_bd_inputs/lsrm_results_v2/vcsn_80perc.h5".format(smt.sdp) # uses the modified 2012 irrgation layer
    outpath = os.path.join(os.path.dirname(path), 'wym_{}'.format(os.path.basename(path)))
    outdata = gen_water_year_average_lsr_irr(path)
    outdata.to_hdf(outpath, 'wym', mode='w')
    rch = map_rch_to_array(hdf=outpath,
                           method='mean',
                           period_center=None,
                           mapping_shp="{}/m_ex_bd_inputs/lsrm_results_v2/test/output_test2.shp".format(smt.sdp),
                           period_length=None, return_irr_demand=False,
                           rch_quanity='total_drainage')

    # fix tewai and chch weirdeness
    fixer = get_rch_fixer()
    # chch
    rch[fixer == 0] = 0.0002
    # te wai and coastal
    rch[fixer == 1] = 0

    # set ibound to 0
    rch[~new_no_flow[0].astype(bool)] = 0
    pickle.dump(rch, open(pickle_path, 'w'))
    return rch


def get_rch_fixer(recalc=False):
    """
    an array to index the abnormal recharge in chch and te waihora 1 = tewaihora and coastal, 0 = chch, all others nan
    :param recalc: boolean whether to recalc (True) or load from pickle if avalible
    :return:
    """
    pickle_path = os.path.join(smt.pickle_dir, 'rch_fixer.p')

    if os.path.exists(pickle_path) and not recalc:
        return pickle.load(open(pickle_path))

    fixer = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/rch_rm_chch_tew.shp".format(smt.sdp), 'ID', True)
    pickle.dump(fixer, open(pickle_path, 'w'))
    return fixer


def _get_rch_comparison():
    """
    create a budget of the recharge between teh current, pc5 and naturalised states
    :return:
    """
    new_no_flow = smt.get_no_flow()

    paths = {'pc5': 'vcsn_100perc.h5', 'cur': 'vcsn_80perc.h5', 'nat': 'vcsn_no_irr.h5'}
    outdict = {}
    for key in paths:
        path = "{}/m_ex_bd_inputs/lsrm_results_v2/{}".format(smt.sdp, paths[key])
        outpath = os.path.join(os.path.dirname(path), 'wym_{}'.format(os.path.basename(path)))
        outdata = gen_water_year_average_lsr_irr(path)
        outdata.to_hdf(outpath, 'wym', mode='w')
        rch = map_rch_to_array(hdf=outpath,
                               method='mean',
                               period_center=None,
                               mapping_shp="{}/m_ex_bd_inputs/lsrm_results_v2/test/output_test2.shp".format(smt.sdp),
                               period_length=None, return_irr_demand=False,
                               rch_quanity='total_drainage')

        # fix tewai and chch weirdeness
        fixer = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/rch_rm_chch_tew.shp".format(smt.sdp), 'ID', True)
        # chch
        rch[fixer == 0] = 0.0002
        # te wai and coastal
        rch[fixer == 1] = 0

        # set ibound to 0
        rch[~new_no_flow[0].astype(bool)] = 0

        outdict[key] = rch

        irr_idx = np.isfinite(
            smt.shape_file_to_model_array("{}\m_ex_bd_inputs\shp\wai_irr_area_intersect.shp".format(smt.sdp), 'year_irr',
                                          'True'))
    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = np.nan
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    w_idx = np.isclose(zones, 4)
    c_idx = np.isclose(zones, 7)
    s_idx = np.isclose(zones, 8)
    all_idx = np.isfinite(zones)
    outdata = pd.DataFrame(columns=['waimak', 'selwyn', 'chch_wm', 'total'])
    for key in outdict:
        dat = outdict[key] * 200 * 200
        for idx, zone in zip([w_idx, c_idx, s_idx, all_idx], ['waimak', 'chch_wm', 'selwyn', 'total']):
            outdata.loc[key, zone] = np.nansum(dat[idx])
    print(outdata / 86400)
    print('done')
    return outdata



if __name__ == '__main__':
    # tests
    test_type = 1


    if test_type == 2:
        rch = _get_rch(2)
        new_no_flow = smt.get_no_flow()
        zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp),'ZONE_CODE')
        zones[~new_no_flow[0].astype(bool)] = np.nan
        # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
        w_idx = np.isclose(zones,4)
        c_idx = np.isclose(zones,7)
        s_idx = np.isclose(zones,8)
        all_idx = np.isfinite(zones)
        irr_idx = np.isfinite(smt.shape_file_to_model_array("{}\m_ex_bd_inputs\shp\wai_irr_area_intersect.shp".format(smt.sdp),'year_irr','True'))
        print('done')


    if test_type ==1:
        rch=_get_rch(version=2,recalc=False)
        x,y = smt.get_model_x_y()
        idx = np.isfinite(rch)
        outdata = pd.DataFrame({'x':x[idx],'y':y[idx],'rch_mm_yr':rch[idx]*1000*365,'rch_m3_s':rch[idx]*200*200/86400})
        outdata.to_csv(os.path.join(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va",'rch_points.csv'))
        raise
        _get_rch_comparison()
        rchold=_get_rch(version=1,recalc=False)

        smt.plt_matrix(rch)
        np.savetxt(r"C:\Users\MattH\Desktop\to_brioch_2017_10_4\rch.txt",rch)
        print('done')