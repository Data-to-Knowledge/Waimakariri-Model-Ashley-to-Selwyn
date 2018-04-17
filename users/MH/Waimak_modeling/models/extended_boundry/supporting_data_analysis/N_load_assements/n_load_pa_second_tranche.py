# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 29/03/2018 3:57 PM
"""

from __future__ import division
from core import env
import geopandas as gpd
import os
from n_load_pa_rules import create_farm_scale_data
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.modified_nload_uncertainty import \
    run_all_nload_stuffs
from sumarize_pa_rules import sumaraize


def run_nload_pa_stuffs(shp_dir, outdir, name, output_act_n=True):
    print('################# {} ################'.format(name))
    run_all_nload_stuffs(base_outdir=os.path.join(outdir, 'stocastic_n_{}'.format(name)), szdirs=[shp_dir],
                         output_act_n=output_act_n)

    # pa stuff
    cments = {}
    shp_paths = glob(os.path.join(shp_dir, "*.shp"))
    for pth in shp_paths:
        cments[os.path.basename(pth).replace('.shp', '')] = gpd.read_file(pth)
    create_farm_scale_data(catchments=cments, outdir=os.path.join(outdir, 'pa_rules_{}'.format(name)))
    sumaraize(os.path.join(outdir, 'pa_rules_{}'.format(name)))


if __name__ == '__main__':
    base_outdir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results"
    base_shp_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons"
    if False:  # quick thing to record, but not re-run
        # private wells
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'private_wells'),
                            outdir=os.path.join(base_outdir, 'private_wells'),
                            name='private_wells')
        # wdc wells
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'wdc_wells'),
                            outdir=os.path.join(base_outdir, 'wdc_wells'),
                            name='wdc_wells')
        # waimakariri
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'waimakariri_river'),
                            outdir=os.path.join(base_outdir, 'waimakariri_river'),
                            name='waimakariri_river')
        # interzone
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'interzone'),
                            outdir=os.path.join(base_outdir, 'interzone'),
                            name='interzone')
        # springfed streams
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'second_tranche'),
                            outdir=os.path.join(base_outdir, 'nwaimak_springfeds'),
                            name='nwaimak_springfeds') #spring fed streams
        # wdc_wells 90 perentage (smaller area)
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'wdc_wells_90_named_right'),
                            outdir=os.path.join(base_outdir, 'wdc_wells_90'),
                            name='wdc_wells_90') #spring fed streams

        # private wells 90 perentage (smaller area)
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'private_wells_90_named_right'),
                            outdir=os.path.join(base_outdir, 'private_wells_90'),
                            name='private_wells_90') #spring fed streams

    if True:
        # wdc wells 90 perentage (smaller area) for all but those over 1000 people
        run_nload_pa_stuffs(shp_dir=os.path.join(base_shp_dir, 'wdc_use_mix'),
                            outdir=os.path.join(base_outdir, 'wdc_use_mix'),
                            name='wdc_use_mix') #spring fed streams
