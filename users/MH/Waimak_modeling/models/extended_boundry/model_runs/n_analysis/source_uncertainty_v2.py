# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 29/03/2018 3:31 PM
"""

from __future__ import division
from core import env
from run_source_uncertainty import calc_all_ns, output_actual_n_vals
import os
import pandas as pd

def run_all_nload_stuffs():
    """
    wrapper to run all the nloads stocastics
    :return:
    """
    base_outdir = env.gw_met_data(r"mh_modeling\stocastic_n_load_results\second_tranche")
    if not os.path.exists(base_outdir):
        os.makedirs(base_outdir)

    szdirs = [
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche"
    ]
    with open(env.gw_met_data(r"mh_modeling\stocastic_n_load_results\n_load_names.txt")) as f:
        n_names = []
        for e in f.readlines():
            if e.strip() != '':
                n_names.append(e.strip())
    for sim_end in ['without_trans', 'with_trans']: # with and without transition to CMP
        for sz_dir in szdirs:
            for n_name in n_names:
                print('starting N analysis for {} load, {} sims, and {} polygons'.format(n_name, sim_end,
                                                                                         os.path.basename(sz_dir)))
                outdir = os.path.join(base_outdir, sim_end, '{}_{}'.format(n_name, os.path.basename(sz_dir)))
                sims = pd.read_csv(env.gw_met_data(
                    "mh_modeling\stocastic_n_load_results\component_uncertainty_data_{}.csv".format(sim_end)),
                    index_col=0)
                calc_all_ns(sims_org=sims, n_load_name=n_name, outdir=outdir, source_zone_dir=sz_dir)
                output_actual_n_vals(outdir=outdir, mod_dir=outdir)

if __name__ == '__main__':
    run_all_nload_stuffs()