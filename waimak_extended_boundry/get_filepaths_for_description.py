# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/07/2018 8:56 AM
"""

from __future__ import division
import env
import os
import pandas as pd
import netCDF4 as nc

# starting with model simulations and results
paths = [
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va",

    # first recursive
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\forward_sw_gw",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\grid_sd",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\lsr_checks",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols",

    # second recursive
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons",

    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\forward_sw_gw\layer0_full_abs_dds",
    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\forward_sw_gw\results",

    r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\lsr_checks\maps",
]

paths2 = [

    r"K:\mh_modeling",

    # first iteration

    r"K:\mh_modeling\data_from_gns",
    r"K:\mh_modeling\data_to_mark",
    r"K:\mh_modeling\data_to_mark.zip",
    r"K:\mh_modeling\interzone_source_zones",
    r"K:\mh_modeling\netcdfs_of_key_modeling_data",
]

out = []
for path in paths2 + paths:
    try:
        out.extend([os.path.join(path, e) for e in os.listdir(path)])
    except:
        pass

# for e in out:
#    print('r"{}",'.format(e))

out = pd.DataFrame({'paths':sorted(out), 'description':[None for e in out], 'Tag':[None for e in out]})

for i in out.index:
    p = out.loc[i,'paths']
    if '.nc' in p:
        temp = nc.Dataset(p)
        out.loc[i,'description'] = temp.description


out.to_csv('c:/users/matth/downloads/paths_to_data.csv', index=False)
