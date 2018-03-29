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

if __name__ == '__main__':
    cments = {}
    shp_paths = glob(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\*.shp")
    for pth in shp_paths:
        cments[os.path.basename(pth).replace('.shp','')] = gpd.read_file(pth)
    create_farm_scale_data(catchments=cments, outdir=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\pa_rules_second_trance")