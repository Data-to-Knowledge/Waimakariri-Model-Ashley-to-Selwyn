# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 21/12/2017 10:52 AM
"""

from __future__ import division
from waimak_extended_boundry.model_run_tools.modpath_source_zone_support.source_raster import get_all_cbcs,run_forward_emulators, get_modeflow_dir_for_source
import os


def run_all_forward_emulators(nsmc_nums, notes, base_results_dir, other_model_ids=None, minparts=1, maxparts=100,
                              create_weak_sink_emulators=True, create_strong_sink_emulators=True):
    """
    run all the forward emulators a big wrapper
    :param nsmc_nums: the nsmc numbers
    :param notes: note to write to readmes
    :param base_results_dir: where to put the results
    :param other_model_ids: other model ids to run (e.g. text model ids)
    :param minparts: the minimum number of particles per cell
    :param maxparts: the maximum number of particles per cell
    :param create_weak_sink_emulators: Bool if True then run for the weak sink emulators
    :param create_strong_sink_emulators: bool if True the run for the strong sink emulators
    :return:
    """
    notes = notes + ' minimum particles: {} maximum particles: {}'.format(minparts, maxparts)
    model_ids = ['NsmcReal{:06d}'.format(e) for e in nsmc_nums]
    if other_model_ids is not None:
        assert isinstance(other_model_ids, list)
        model_ids.extend(other_model_ids)
    modflow_dir = get_modeflow_dir_for_source()
    strong_results_dir = os.path.join(base_results_dir, 'strong_sinks')
    weak_results_dir = os.path.join(base_results_dir, 'weak_sinks')

    # check all cbcs have been created
    get_all_cbcs(model_ids, modflow_dir)
    print('starting to run particle tracking')
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
    raise ValueError('I probably dont want to re-run this')
    #nsmc_nums = [5, 17, 18, 26, 37, 44, 72, 103, 117, 133, 142,
    #             204, 233, 240, 258, 271, 278, 308, 314, 388, 391, 439,
    #             441, 447, 488, 491, 552, 555, 604, 640, 666, 790, 794,
    #             800, 808, 809, 871, 883, 932, 952, 955, 992, 1001, 1029,
     #            1052, 1055, 1078, 1099, 1106, 1111, 1112, 1114, 1129, 1168, 1181,
      #           1250, 1326, 1328, 1330, 1349, 1365, 1368, 1370, 1395, 1399, 1454,
       #          1481, 1489, 1498, 1499, 1507, 1574, 1585, 1595, 1618, 1636, 1656,
        #         1660, 1681, 1694, 1749, 1773, 1817, 1819, 1836, 1842, 1851, 1870,
         #        1881, 1897, 1915, 1916, 1967, 2044, 2052, 2068, 2085, 2114, 2117,
          #       2139, 2160, 2200, 2209, 2234, 2278, 2382, 2390, 2430, 2460, 2520,
           #      2579, 2588, 2590, 2604, 2643, 2654, 2661, 2673, 2691, 2735, 2748,
            #     2752, 2757, 2760, 2818, 2840, 2871, 2875, 2885, 2918, 2971, 2980,
             #    2997, 3039, 3064, 3067, 3075, 3100, 3116, 3165, 3216, 3289, 3310,
              #   3318, 3350, 3378, 3389, 3418, 3428, 3456, 3464, 3482, 3513, 3585,
               #  3641, 3651, 3663, 3685, 3712, 3717, 3762, 3796, 3836, 3861, 3910] # the 165 models which passed the EMMA filter
    nsmc_nums = [2971, 3116, 3310, 3378, 3456, 3513, 3762, 3910] # rerunning errrors
    notes = 'forward runs for the 165 models plus the ashley river optimisation for the calibration period'
    base_results_dir = get_base_results_dir('forward', socket.gethostname())
    #other_model_ids = ['NsmcBase', 'AshOpt', 'StrOpt']
    other_model_ids = [] # re-running some that missed
    run_all_forward_emulators(nsmc_nums, notes, base_results_dir,
                              other_model_ids=other_model_ids,
                              minparts=1, maxparts=20,
                              create_strong_sink_emulators=True,
                              create_weak_sink_emulators=True)
