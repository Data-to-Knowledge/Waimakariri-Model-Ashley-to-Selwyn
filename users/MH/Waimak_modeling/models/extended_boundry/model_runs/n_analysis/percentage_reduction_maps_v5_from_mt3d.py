"""
Author: matth
Date Created: 3/05/2018 1:34 PM
"""

from __future__ import division
from core import env
import itertools
from percentage_reduction_maps import gen_stream_targets, gen_well_targets,gen_waimak_targets,\
    calc_per_reduction_rasters, get_mode, get_interzone_reduction


if __name__ == '__main__':
    #todo set up new scenario, and setup from mt3d run options
    #todo what about the interzone reductions, do we want them?
    test = get_interzone_reduction(8,True)
    test2 = get_interzone_reduction(8,False)
    outdir = (r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_res"
             r"ults\n_reductions_from_interzone")
    interzone_target_load = 8 #todo I can also pass a '50%' str

    # 3 main scenarios
    scenarios = ['least_pain', 'middle_option', 'most_gain'] #todo note that the mode is defined for each scenario
    # with and without and pc5pa
    mar_pc5pa = [True, False]
    # with mode  = 50th and 95th # todo once from mt3d I should be able to pass modes of any variety (e.g. standard naming)
    # with and without conservative things
    conservative_shps = ['use_mix']

    for scen, mar, conserv in itertools.product(scenarios,mar_pc5pa, conservative_shps):
        print(scen, mar, conserv)
        mode = get_mode(scen, current_from_mt3d=True)

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
                                   interzone_target_load=interzone_target_load, include_interzone=True, current_pathway_from_mt3d=True)

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
                                   interzone_target_load=None, include_interzone=False, save_reason=False, current_pathway_from_mt3d=True)

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
                                   interzone_target_load=None, include_interzone=False, save_reason=False, current_pathway_from_mt3d=True)

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
                                   interzone_target_load=None, include_interzone=False, save_reason=False, current_pathway_from_mt3d=True)

        # interzone only
        well_targets = gen_well_targets(scen,wdc_none=True, private_none=True)
        stream_targets = gen_stream_targets(scen, True)
        waimak_target = 0
        if mar:
            mar_per = 0  # assume no mar
        else:
            mar_per = 0
        calc_per_reduction_rasters(outdir, '{}_interzone_only'.format(name), mode, well_targets, stream_targets, waimak_target=waimak_target,
                                   mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                   interzone_target_load=interzone_target_load, include_interzone=True, save_reason=False, current_pathway_from_mt3d=True)

