# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 21/12/2017 10:52 AM
"""

from __future__ import division
from core import env
from source_raster import get_all_cbcs,run_forward_emulators, get_modeflow_dir_for_source
import os

if __name__ == '__main__':
    # items to change for different runs
    nsmc_nums = [-1,-2] #todo change
    notes = 'first try' #todo change
    base_results_dir = r"D:\mh_waimak_models\modpath_test_forward" #todo change make sure this is on the C drive as it will be much faster

    minparts = 1
    maxparts = 20


    create_cbcs = False
    create_weak_sink_emulators = False
    create_strong_sink_emulators = True
    # these items end here

    model_ids = ['NsmcReal{:06d}'.format(e) for e in nsmc_nums]
    modflow_dir = get_modeflow_dir_for_source()
    strong_results_dir = os.path.join(base_results_dir,'strong_sinks')
    weak_results_dir = os.path.join(base_results_dir,'weak_sinks')

    if create_cbcs:
        get_all_cbcs(model_ids,modflow_dir)
    if create_weak_sink_emulators:
        # 30 minutes for min part 1 max part 100
        # 45 minutes for min part 1 max part 500
        run_forward_emulators(model_ids, weak_results_dir, modflow_dir, keep_org_files=False, min_part=minparts,
                              max_part=maxparts, capt_weak_s=True, notes=notes)

    if create_strong_sink_emulators:
        run_forward_emulators(model_ids, strong_results_dir, modflow_dir, keep_org_files=False, min_part=minparts,
                              max_part=maxparts, capt_weak_s=False, notes=notes)
    #todo debug
