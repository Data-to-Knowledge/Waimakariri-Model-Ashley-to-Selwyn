# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 21/12/2017 10:52 AM
"""

from __future__ import division
from core import env
from source_raster import get_all_cbcs,run_forward_emulators, get_modeflow_dir_for_source
import os

def run_all_forward_emulators(nsmc_nums, notes, base_results_dir, other_model_ids=None, minparts=1, maxparts=100,
                              create_weak_sink_emulators=True, create_strong_sink_emulators=True):
    notes = notes + ' minimum particles: {} maximum particles: {}'.format(minparts,maxparts)
    model_ids = ['NsmcReal{:06d}'.format(e) for e in nsmc_nums]
    if other_model_ids is not None:
        assert isinstance(other_model_ids,list)
        model_ids.extend(other_model_ids)
    modflow_dir = get_modeflow_dir_for_source()
    strong_results_dir = os.path.join(base_results_dir,'strong_sinks')
    weak_results_dir = os.path.join(base_results_dir,'weak_sinks')

    # check all cbcs have been created
    get_all_cbcs(model_ids,modflow_dir)
    if create_weak_sink_emulators:
        # 30 minutes for min part 1 max part 100
        # 45 minutes for min part 1 max part 500
        run_forward_emulators(model_ids, weak_results_dir, modflow_dir, keep_org_files=False, min_part=minparts,
                              max_part=maxparts, capt_weak_s=True, notes=notes)

    if create_strong_sink_emulators:
        run_forward_emulators(model_ids, strong_results_dir, modflow_dir, keep_org_files=False, min_part=minparts,
                              max_part=maxparts, capt_weak_s=False, notes=notes)


if __name__ == '__main__':
    # items to change for different runs
    nsmc_nums = [-1,-2]  # todo change
    notes = 'forward runs for the 165 models plus the ashley river optimisation for the calibration period'
    base_results_dir = r"C:\mh_waimak_models\modpath_forward_base"
    other_model_ids = None  # could add a list of other model ids
    run_all_forward_emulators(nsmc_nums,notes,base_results_dir,
                              other_model_ids=other_model_ids,
                              minparts=1, maxparts=20,
                              create_strong_sink_emulators=True,
                              create_weak_sink_emulators=True)
