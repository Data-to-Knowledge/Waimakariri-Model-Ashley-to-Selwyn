# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/10/2017 4:14 PM
"""

from __future__ import division
import os
from waimak_extended_boundry import \
    get_full_consent, get_max_rate
from base_modflow_wrapper import import_gns_model
from realisation_id import get_base_well, get_rch_multipler, get_model_name_path
from waimak_extended_boundry import \
    get_base_grid_sd_path
from waimak_extended_boundry import \
    get_str_dep_base_path
from waimak_extended_boundry import \
    _get_no_pumping_ss_hds
from waimak_extended_boundry import modflow_converged
from waimak_extended_boundry import \
    get_drn_no_ncarpet_spd

def minium_new_model_run(model_id):
    test_path = get_model_name_path(model_id)
    get_base_well(model_id, org_pumping_wells=True, recalc=True)
    get_base_well(model_id, org_pumping_wells=False, recalc=True)
    get_full_consent(model_id, org_pumping_wells=True, recalc=True)
    get_full_consent(model_id, org_pumping_wells=False, recalc=True)
    get_max_rate(model_id, org_pumping_wells=True, recalc=True)
    get_max_rate(model_id, org_pumping_wells=False, recalc=True)
    get_drn_no_ncarpet_spd(model_id, recalc=True)

    print('can start running models')



if __name__ == '__main__':
    # a script to pull through all of the pickles needed for a new model_id
    model_id = 'NsmcReal{:06d}'.format(491)


    t = False  # do not change instead change indentation
    if t:  # just a quick way to skip things when debugging
        print('stuff to skip')
        get_model_name_path(model_id)
        get_rch_multipler(model_id)
        m = import_gns_model(model_id, 'test', r"C:\Users\MattH\Desktop\test")
        print(m)
        m.write_input()
        m.write_name_file()
        m.run_model()
        print('{} at least loaded'.format(model_id))
        con = modflow_converged(os.path.join(m.model_ws, '{}.list'.format(m.name)))
        if not con:
            raise ValueError('model did not converge')
        get_str_dep_base_path(model_id, 'sd7', recalc=True)
        get_str_dep_base_path(model_id, 'sd30', recalc=True)
        get_str_dep_base_path(model_id, 'sd150', recalc=True)
        get_base_grid_sd_path(model_id, recalc=True)
        _get_no_pumping_ss_hds(model_id, recalc=True)

    get_full_consent(model_id, org_pumping_wells=True, recalc=True)
    get_full_consent(model_id, org_pumping_wells=False, recalc=True)
    get_max_rate(model_id, org_pumping_wells=True, recalc=True)
    get_max_rate(model_id, org_pumping_wells=False, recalc=True)
    get_base_well(model_id, org_pumping_wells=True, recalc=True)
    get_base_well(model_id, org_pumping_wells=False, recalc=True)

    print('can start running models')
    get_drn_no_ncarpet_spd(model_id, recalc=True)

