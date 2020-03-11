# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 26/10/2017 3:22 PM
"""
#todo it might make sense to move this into the model run tools, it is used a lot
from __future__ import division
from waimak_extended_boundry.model_and_NSMC_build.m_packages.drn_packages import _get_drn_spd
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools import \
    mod_gns_model, get_model, get_race_data
import numpy as np
import os
import pickle
import flopy_mh as flopy

def get_drn_spd(): #todo move this up
    raise NotImplementedError

def get_drn_no_ncarpet_spd(model_id,recalc=False):
    pickle_path = "{}/model_{}_drn_wout_ncarpet.p".format(smt.temp_pickle_dir, model_id)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = pickle.load(open(pickle_path))
        return outdata
    model = get_model(model_id)
    org_drn_data = _get_drn_spd(1, 1)
    n_carpet = org_drn_data.loc[np.in1d(org_drn_data.group, ['ash_carpet', 'cust_carpet'])]
    n_carpet_ij = list(n_carpet.loc[:, ['i', 'j']].itertuples(False, None))
    outdata = model.get_package('drn').stress_period_data.data[0]
    idx = [e not in n_carpet_ij for e in zip(outdata['i'], outdata['j'])]

    outdata = outdata[idx]
    pickle.dump(outdata,open(pickle_path,'w'))

    return outdata

def raise_top_elv(model):
    # try first only in carpet drain cells, then second elsewhere
    raise NotImplementedError


if __name__ == '__main__':
    # look for difference with and without carpet drains and/or raised elevations
    test_type = [4]
    if 1 in test_type:
        model_id = 'StrOpt'
        m = mod_gns_model(model_id, 'test_raise_elv',
                          r"C:\Users\MattH\Desktop\test_raise_elv_raised", False,drain={0:get_drn_no_ncarpet_spd(model_id)})
        m.dis.top= flopy.utils.Util2d(m,m.dis.top.shape,m.dis.top.dtype,m.dis.top.array+10,m.dis.top.name)
        m.write_input()
        m.write_name_file()
        m.run_model()
    if 2 in test_type:
        hds = flopy.utils.HeadFile(r"C:\Users\MattH\Desktop\StrOpt_test_raise_elv\StrOpt_test_raise_elv.hds").get_data(kstpkper=(0,0))
        hds[hds>1e10] = np.nan
        elv = smt.calc_elv_db()
        temp = hds[0]-elv[0]

    if 3 in test_type:
        hds_cal = flopy.utils.HeadFile(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\StrOpt\AW20171019_i3_optver\mf_aw_ex.hds").get_data(kstpkper=(0,0))
        hds_cal[hds_cal>1e10] = np.nan

        hds_org = flopy.utils.HeadFile(r"C:\Users\MattH\Desktop\StrOpt_test_raise_elv\StrOpt_test_raise_elv.hds").get_data(kstpkper=(0,0))
        hds_org[hds_org>1e10] = np.nan

        hds_raised = flopy.utils.HeadFile(r"C:\Users\MattH\Desktop\StrOpt_test_raise_elv_raised\StrOpt_test_raise_elv.hds").get_data(kstpkper=(0,0))
        hds_raised[hds_raised>1e10] = np.nan
    if 4 in test_type:
        model_id = 'StrOpt'
        wells = {}
        base_well = get_race_data(model_id)
        wells[0] = smt.convert_well_data_to_stresspd(base_well)

        m = mod_gns_model(model_id, 'test_raise_elv',
                          r"C:\Users\MattH\Desktop\test_raise_elv_raised_nat", False,
                          drain={0:get_drn_no_ncarpet_spd(model_id)},
                          well=wells)
        m.dis.top= flopy.utils.Util2d(m,m.dis.top.shape,m.dis.top.dtype,m.dis.top.array+10,m.dis.top.name)
        m.write_input()
        m.write_name_file()
        m.run_model()




    print('done')
