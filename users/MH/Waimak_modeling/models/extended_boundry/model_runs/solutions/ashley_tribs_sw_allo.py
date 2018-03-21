# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 21/03/2018 9:56 AM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.base_modflow_wrapper import \
    mod_gns_model
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import get_wel_spd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.stream_depletion_assesment.raising_heads_no_carpet import \
    get_drn_no_ncarpet_spd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.min_flows_reliability.extract_data_for_forward_runs import \
    extract_and_save_all_forward_runs
import multiprocessing
import pandas as pd
import geopandas as gpd
import os
import psutil
import logging
import numpy as np

def setup_run_new_abstraction(name, new_wells, base_path):
    drn_data = get_drn_no_ncarpet_spd('NsmcBase')
    well_data = get_wel_spd(3)
    if new_wells is not None:
        well_data = pd.concat((well_data, new_wells))
    m = mod_gns_model('NsmcBase', name, os.path.join(base_path, name),
                      False, well={0:smt.convert_well_data_to_stresspd(well_data)},
                      drain={0:drn_data})
    m.write_input()
    m.run_model()


def setup_run_new_abstraction_mp(kwargs):
    setup_run_new_abstraction(**kwargs)


def start_process():
    """
    function to run at the start of each multiprocess sets the priority lower
    :return:
    """
    print('Starting', multiprocessing.current_process().name)
    p = psutil.Process(os.getpid())
    # set to lowest priority, this is windows only, on Unix use ps.nice(19)
    p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)


def get_new_wells():
    new = gpd.read_file(r'{}\m_ex_bd_inputs\shp\solutions_ashley_new_wells.shp'.format(smt.sdp))
    xs = new.geometry.x
    ys = new.geometry.y
    rows = np.nan * xs
    cols = np.nan * ys
    # put in layer 5
    layer = ys * 0 + 5
    flux = ys * 0 + 10 * 86.4  # 10 l/s converted to m to day
    for i in range(len(xs)):
        r, c = smt.convert_coords_to_matix(xs[i], ys[i])
        rows[i] = r
        cols[i] = c

    outdata = pd.DataFrame({'stream': new['stream'], 'col': cols, 'row': rows, 'layer': layer, 'flux': flux})
    return outdata

def main():
    all_new_well = get_new_wells()
    base_dir = r"C:\Users\MattH\Downloads\sw_solutions_ashley"
    outpath = os.path.join(base_dir,'solutions.csv')
    readme_txt = """the following runs are done to examine the potential of moving some surface water abstraction 
    into deep wells.  10l deep wells have been added to the following swaps saltwater(2 well) taranaki(9 wells) 
    waikuku(14 wells) as indicated by names (e.g. saltwater only has saltwater wells added). 
    Carpet drains have been removed""".replace('\n', '')
    runs = []
    runs.append({'name': 'original', 'new_wells': None, 'base_path': base_dir})
    runs.append({'name': 'all', 'new_wells': all_new_well, 'base_path': base_dir})
    # make kwargs
    for stream in set(all_new_well.stream):
        runs.append({'name': stream,
                     'new_wells': all_new_well.loc[all_new_well.stream==stream],
                     'base_path': base_dir})


    multiprocessing.log_to_stderr(logging.DEBUG)
    pool_size = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(processes=pool_size,
                                initializer=start_process,
                                )
    results = pool.map(setup_run_new_abstraction_mp, runs)
    extract_and_save_all_forward_runs(base_dir, outpath, readme_txt=readme_txt)
    with open(os.path.join(base_dir, 'READ_ME.txt'), 'w') as f:
        f.write(readme_txt)


if __name__ == '__main__':
    main()
