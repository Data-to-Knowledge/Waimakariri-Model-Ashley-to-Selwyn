# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/10/2017 4:48 PM
"""

from __future__ import division
from core import env
import os
import time
from well_by_well_str_dep_sd7 import get_sd_well_list, well_by_well_depletion_sd7
from well_by_well_str_dep_sd30 import well_by_well_depletion_sd30
from well_by_well_str_dep_sd150 import well_by_well_depletion_sd150
from extract_stream_depletion import calc_str_dep_all_wells
import shutil
from sd_metadata import save_sd_metadata
from starting_hds_ss_sy import get_ss_sy
import datetime
import socket


def run_full_non_grid_sd(model_id, model_base_path, data_outdir, ss_sy_version, notes=''):
    """
    a conveninance function to run all of the non-well stream depletion assesments and data extraction Note that this
    does not keep the model file, those are deleted after data extraction
    :param model_id: which nsmc realisation to use
    :param model_base_path: the base path to put the models
    :param data_outdir: the base path to put the data
    :param ss_sy_version: which version of ss_sy to use see get_ss_sy()
    :param notes: any notes to write to a readme file
    :return:
    """

    description = get_ss_sy(ss_sy_version=ss_sy_version, return_description=True)
    sd7_notes = sd30_notes = sd150_notes = notes + '\nss_sy version: {}\n ss_sy_description:{}'.format(ss_sy_version,
                                                                                                       description)
    # below should not change
    if not os.path.exists(model_base_path):
        os.makedirs(model_base_path)
    if not os.path.exists(data_outdir):
        os.makedirs(data_outdir)

    with open(os.path.join(data_outdir, 'READ_ME.txt'), 'w') as f:
        f.write(sd7_notes)

    well_list = ['M35/10410','M35/0335', 'M35/4238'] #get_sd_well_list(model_id) # note I'm re-running some wells that missed because of bug

    #### run the models, extract data, delete model files ####
    sd7_base_path = os.path.join(model_base_path, 'sd7')
    sd30_base_path = os.path.join(model_base_path, 'sd30')
    sd150_base_path = os.path.join(model_base_path, 'sd150')

    t = time.time()
    # sd 7
    print('starting to run SD7 models')  # rdsprod 15 min gw02: ~ 25 min for 10
    well_by_well_depletion_sd7(model_id, well_list, sd7_base_path, sd7_notes, ss_sy_version=ss_sy_version)

    print('extracting data from SD7')
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd7.csv"), sd7_base_path, 'sd7', ss_sy_version)
    save_sd_metadata(os.path.join(data_outdir, "sd7_metadata.csv"), sd7_base_path)

    print('removing model files for SD7')
    shutil.rmtree(sd7_base_path)

    # sd 30
    print('starting to run SD30 models')  # rdsprod 15 min gw02 ~ 35 min for 10 models
    well_by_well_depletion_sd30(model_id, well_list, sd30_base_path, sd30_notes, ss_sy_version=ss_sy_version)

    print('extracting data from SD30')
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd30.csv"), sd30_base_path, 'sd30', ss_sy_version)
    save_sd_metadata(os.path.join(data_outdir, "sd30_metadata.csv"), sd30_base_path)

    print('removing model files for SD30')
    shutil.rmtree(sd30_base_path)

    # sd 150
    print('starting to run SD150 models')  # rdsprod 45 gw02 ~ 100 min for 10 models
    well_by_well_depletion_sd150(model_id, well_list, sd150_base_path, sd150_notes, ss_sy_version=ss_sy_version)

    print('extracting data from SD150')
    calc_str_dep_all_wells(os.path.join(data_outdir, "extract_sd150.csv"), sd150_base_path, 'sd150', ss_sy_version)
    save_sd_metadata(os.path.join(data_outdir, "sd150_metadata.csv"), sd150_base_path)

    print('removing model files for SD150')
    shutil.rmtree(sd150_base_path)

    print('all sd  done after {} minutes for {} model runs'.format((time.time() - t) / 60, len(well_list * 3)))
    print(notes)


if __name__ == '__main__':
    # well list is 732 models
    # ~ 5 days per storage value for rdsprod 03 (8 cores)
    # ~ 6.75 days per storage value for gw02 (12 cores)
    # above based on lowest storage value, which is probably the hardest
    model_id = 'NsmcBase'
    if socket.gethostname() == 'RDSProd03':
        sd_sy_versions = [2]
        sd_sy_v_names = ['low_s']
    elif socket.gethostname() == 'GWATER02':
        sd_sy_versions = [3] #4] # see below
        sd_sy_v_names = ['med_s'] # 'high_s'] # high was also run on gw02, but they also killed RDS prod03
        sd_sy_versions = [2,3,4] #re-run a couple missed
        sd_sy_v_names = ['low_s','med_s', 'high_s'] #re-run a couple missed

    else:
        raise ValueError('not set up for {}'.format(socket.gethostname()))

    base_dir = "D:/mh_waimak_models/sd_multi_s_re-run"
    base_notes = 'a test of a couple of wells through the for loop'

    today = str(datetime.date.today()).replace('-', '_')
    for version, name in zip(sd_sy_versions, sd_sy_v_names):
        version_dir = os.path.join(base_dir, name)

        if not os.path.exists(version_dir):
            os.makedirs(version_dir)

        model_base_path = os.path.join(version_dir, "{}_{}_models".format(model_id, today))
        data_outdir = os.path.join(version_dir, "{}_{}_data".format(model_id, today))

        run_full_non_grid_sd(model_id, model_base_path, data_outdir, ss_sy_version=version,
                             notes=base_notes)
