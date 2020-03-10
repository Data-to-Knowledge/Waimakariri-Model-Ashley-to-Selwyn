# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 16/05/2018 3:13 PM
"""

from __future__ import division
import itertools
from percentage_reduction_maps import gen_stream_targets, gen_well_targets, gen_waimak_targets, \
    calc_per_reduction_rasters, get_mode
from plot_reduction_rasters import plot_all_rasters
import time
import os
from amount_landuse_change import get_percentages
import pandas as pd

outdir = (r"C:\Users\MattH\Downloads\test_in_meeting_reduction_rasters_v2")
scenarios = ['most_gain']  # note that the mode is defined for each scenario

if __name__ == '__main__':
    t = time.time()
    # to add a new scenario I must update scenarios in get well targets, get stream targets, get waimak target and get mode
    # to plot run plot_reduction_rasters
    # to get area run amount_landuse_change
    #todo the mode and setup for stream and well targets could be improved...
    make_components = False
    include_interzone = False
    if include_interzone:
        interzone_target_load = 8  # I can also pass a '50%' str
    else:
        interzone_target_load = None

    # 3 main scenarios
    # with and without and pc5pa
    mar_pc5pa = [True, False]
    # with mode  = 50th and 95th # from mt3d I should be able to pass modes of any variety (e.g. standard naming)
    # with and without conservative things
    conservative_shps = ['use_mix']

    for scen, mar, conserv in itertools.product(scenarios, mar_pc5pa, conservative_shps):
        print(scen, mar, conserv)
        mode = get_mode(scen, current_from_mt3d=True)

        con_name = conserv
        if mar:
            mar_name = 'with_pc5pa00'
        else:
            mar_name = 'without_pc5pa00'

        name = '{}_{}_{}'.format(scen, con_name, mar_name)
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
                                   interzone_target_load=interzone_target_load,
                                   include_interzone=include_interzone, current_pathway_from_mt3d=True)

        if make_components:
            # streams only receptors
            well_targets = gen_well_targets(scen, True, True)
            stream_targets = gen_stream_targets(scen)
            waimak_target = gen_waimak_targets(scen)
            if mar:
                mar_per = 0  # assume no mar
            else:
                mar_per = 0
            calc_per_reduction_rasters(outdir, '{}_stream_only'.format(name), mode, well_targets, stream_targets,
                                       waimak_target=waimak_target,
                                       mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                       interzone_target_load=None, include_interzone=False, save_reason=False,
                                       current_pathway_from_mt3d=True)

            # private wells only
            well_targets = gen_well_targets(scen, wdc_none=True, private_none=False)
            stream_targets = gen_stream_targets(scen, True)
            waimak_target = 0
            if mar:
                mar_per = 0  # assume no mar
            else:
                mar_per = 0
            calc_per_reduction_rasters(outdir, '{}_private_wells_only'.format(name), mode, well_targets, stream_targets,
                                       waimak_target=waimak_target,
                                       mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                       interzone_target_load=None, include_interzone=False, save_reason=False,
                                       current_pathway_from_mt3d=True)

            # wdc wells only
            well_targets = gen_well_targets(scen, wdc_none=False, private_none=True)
            stream_targets = gen_stream_targets(scen, True)
            waimak_target = 0
            if mar:
                mar_per = 0  # assume no mar
            else:
                mar_per = 0
            calc_per_reduction_rasters(outdir, '{}_wdc_wells_only'.format(name), mode, well_targets, stream_targets,
                                       waimak_target=waimak_target,
                                       mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                       interzone_target_load=None, include_interzone=False, save_reason=False,
                                       current_pathway_from_mt3d=True)

            # interzone only
            if include_interzone:
                well_targets = gen_well_targets(scen, wdc_none=True, private_none=True)
                stream_targets = gen_stream_targets(scen, True)
                waimak_target = 0
                if mar:
                    mar_per = 0  # assume no mar
                else:
                    mar_per = 0
                calc_per_reduction_rasters(outdir, '{}_interzone_only'.format(name), mode, well_targets, stream_targets,
                                           waimak_target=waimak_target,
                                           mar_percentage=mar_per, pc5_pa_rules=mar, conservative_shp=conserv,
                                           interzone_target_load=interzone_target_load, include_interzone=True,
                                           save_reason=False, current_pathway_from_mt3d=True)

    # plot these things
    print('generating plots')
    plot_all_rasters(outdir,os.path.join(outdir,'plots'))

    # get land reductions
    print('getting land reduction estimates')
    mar_names = ['with_pc5pa00', 'without_pc5pa00']
    # with and without and pc5pa
    # with mode  = 50th and 95th
    # with and without conservative things
    mar_names = []
    for mar in mar_pc5pa:
        if mar:
            mar_names.append('with_pc5pa00')
        else:
            mar_names.append('without_pc5pa00')

    conservative_shps = ['use_mix']

    # mixed
    if make_components:
        appendixes = ['', '_stream_only', '_wdc_wells_only', '_private_wells_only']
    else:
        appendixes = ['']
    for ider in appendixes:
        outdata = pd.DataFrame(index=pd.MultiIndex.from_product([scenarios, mar_names], names=['option', 'scenario']),
                               columns=['Land areas requiring 3-7% N load reduction',
                                        'Land areas requiring 7-20% N load reduction',
                                        'Land areas requiring >20% N load reduction'])
        for scen, mar_name in itertools.product(scenarios, mar_names):
            tif_path = os.path.join(outdir, '{}_use_mix_{}{}_reduction.tif'.format(scen, mar_name, ider))
            print(tif_path)
            temp = get_percentages(tif_path)
            outdata.loc[scen, mar_name] = temp
        outdata *= 100
        outdata = outdata.astype(int)
        outdata.to_csv(os.path.join(outdir, 'areas_mixed{}.csv'.format(ider)))

    print('took {} min'.format((time.time()-t)/60))
