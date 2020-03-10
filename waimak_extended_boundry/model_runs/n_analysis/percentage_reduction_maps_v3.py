# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/04/2018 11:17 AM
"""

from __future__ import division
import itertools
from percentage_reduction_maps import gen_stream_targets, gen_well_targets,gen_waimak_targets,\
    calc_per_reduction_rasters, get_mode


if __name__ == '__main__':
    outdir = (r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_res"
             r"ults\n_reductions_use_zones_use_scens")

    # 3 main scenarios
    scenarios = ['least_pain', 'middle_option', 'most_gain']
    # with and without and pc5pa
    mar_pc5pa = [True, False]
    # with mode  = 50th and 95th
    # with and without conservative things
    conservative_shps = ['use_mix']

    for scen, mar, conserv in itertools.product(scenarios,mar_pc5pa, conservative_shps):
        print(scen, mar, conserv)
        mode = get_mode(scen)

        con_name = conserv
        if mar:
            mar_name = 'with_pc5pa00'
        else:
            mar_name = 'without_pc5pa00'

        name = '{}_{}_{}'.format(scen,con_name,mar_name)
        # all receptors
        well_targets = gen_well_targets(scen)
        stream_targets = gen_stream_targets(scen)
        waimak_target = gen_waimak_targets(scen)
        if mar:
            mar_per = 0  # assume no mar
        else:
            mar_per = 0
        calc_per_reduction_rasters(outdir, name, mode, well_targets, stream_targets, waimak_target=waimak_target,
                                   mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                   interzone_target_load=None, include_interzone=False)

        # streams only receptors
        well_targets = gen_well_targets(scen, True, True)
        stream_targets = gen_stream_targets(scen)
        waimak_target = gen_waimak_targets(scen)
        if mar:
            mar_per = 0  # assume no mar
        else:
            mar_per = 0
        calc_per_reduction_rasters(outdir, '{}_stream_only'.format(name), mode, well_targets, stream_targets, waimak_target=waimak_target,
                                   mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                   interzone_target_load=None, include_interzone=False, save_reason=False)

        # private wells only
        well_targets = gen_well_targets(scen,wdc_none=True, private_none=False)
        stream_targets = gen_stream_targets(scen, True)
        waimak_target = 0
        if mar:
            mar_per = 0  # assume no mar
        else:
            mar_per = 0
        calc_per_reduction_rasters(outdir, '{}_private_wells_only'.format(name), mode, well_targets, stream_targets, waimak_target=waimak_target,
                                   mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                   interzone_target_load=None, include_interzone=False, save_reason=False)

        # wdc wells only
        well_targets = gen_well_targets(scen,wdc_none=False, private_none=True)
        stream_targets = gen_stream_targets(scen, True)
        waimak_target = 0
        if mar:
            mar_per = 0  # assume no mar
        else:
            mar_per = 0
        calc_per_reduction_rasters(outdir, '{}_wdc_wells_only'.format(name), mode, well_targets, stream_targets, waimak_target=waimak_target,
                                   mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                   interzone_target_load=None, include_interzone=False, save_reason=False)

