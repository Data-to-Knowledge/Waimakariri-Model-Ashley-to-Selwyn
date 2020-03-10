# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/10/2017 4:48 PM
"""

from __future__ import division
import os
import time
from well_by_well_str_dep_sd7 import get_sd_well_list, well_by_well_depletion_sd7
from well_by_well_str_dep_sd30 import well_by_well_depletion_sd30
from well_by_well_str_dep_sd150 import well_by_well_depletion_sd150
from extract_stream_depletion import calc_str_dep_all_wells
from sd_metadata import save_sd_metadata

if __name__ == '__main__':
    # a conveninance function to run all of the non-well stream depletion assesments and data extraction

    #### update below parameters ####
    # on DHI-runs02
    model_id = 'NsmcBase'

    sd7_notes = """ without carpet drains in the north"""
    sd30_notes = """ without carpet drains in the north"""
    sd150_notes = """ without carpet drains in the north """


    model_base_path = r"E:\mattH\python_wm_runs\sd_runs\{}_2017_10_26".format(model_id)
    data_outdir = r"E:\mattH\python_wm_runs\sd_runs\data_{}_2017_10_26".format(model_id)
    run_models = True

    # below should not change
    if not os.path.exists(model_base_path):
        os.makedirs(model_base_path)
    if not os.path.exists(data_outdir):
        os.makedirs(data_outdir)

    well_list = get_sd_well_list(model_id)

    #### run the models ####
    sd7_base_path = os.path.join(model_base_path, 'sd7')
    sd30_base_path = os.path.join(model_base_path, 'sd30')
    sd150_base_path = os.path.join(model_base_path, 'sd150')

    if run_models:
        t = time.time()

        well_by_well_depletion_sd7(model_id, well_list, sd7_base_path, sd7_notes)
        well_by_well_depletion_sd30(model_id, well_list, sd30_base_path, sd30_notes)
        well_by_well_depletion_sd150(model_id, well_list, sd150_base_path, sd150_notes)

        print('done after {} minutes for {} model runs'.format((time.time() - t) / 60, len(well_list) * 3))

    #### extract data ####
    if not os.path.exists(data_outdir):
        os.makedirs(data_outdir)
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd150.csv"), sd150_base_path, 'sd150')
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd30.csv"), sd30_base_path, 'sd30')
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd7.csv"), sd7_base_path, 'sd7')

    # metadata
    save_sd_metadata(os.path.join(data_outdir, "sd150_metadata.csv"), sd150_base_path)
    save_sd_metadata(os.path.join(data_outdir, "sd30_metadata.csv"), sd30_base_path)
    save_sd_metadata(os.path.join(data_outdir, "sd7_metadata.csv"), sd7_base_path)

