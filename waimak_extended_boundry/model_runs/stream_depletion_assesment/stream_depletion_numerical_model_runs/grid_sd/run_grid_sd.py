# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/10/2017 5:04 PM
"""

from __future__ import division
import os
import time
from well_by_well_str_dep_ss import well_by_well_depletion_grid
from extract_grid_sd import calc_str_dep_all_wells_grid
from visualise_grid_sd import krig_plot_sd_grid
import socket
from waimak_extended_boundry.model_runs.stream_depletion_assesment.stream_depletion_numerical_model_runs.sd_metadata import save_sd_metadata
from make_grid_sd_zones import make_sd_zone
if __name__ == '__main__':
    # a convenience function to run all of the grid sd, extract, and visualise the data
    # size requirements: one run (one well) is ~ 47 MB there are one full grid run is ~85 GB
    # run time requirments:

    #### update the below parameters ####
    # run on GWruns02
    model_id = 'NsmcBase'  # re-define
    if socket.gethostname() == 'RDSProd03':
        base_dir = "D:\mh_model_runs\grid_sd_runs\{}_models_2017_10_27".format(model_id)
        data_out_dir = "D:\mh_model_runs\grid_sd_runs\{}_data_2017_10_27".format(model_id)
        run_models = True
        amalg_results = True

        # below should not change
        fluxes = [-100, - 25]
        fluxes = [e * 86.4 for e in fluxes]

    else: # runs on gw02
        run_models = True
        amalg_results = True
        base_dir = "D:\mh_waimak_models\grid_sd_runs\{}_models_2017_10_28".format(model_id)
        data_out_dir = "D:\mh_waimak_models\grid_sd_runs\{}_data_2017_10_28".format(model_id)

        # below should not change
        fluxes = [-25]
        fluxes = [e * 86.4 for e in fluxes]

    #### run the models ####
    if run_models:
        t = time.time()
        base_paths = []
        for flux in fluxes:
            path = os.path.join(base_dir, 'flux_{}'.format(flux))
            base_paths.append(path)
            well_by_well_depletion_grid(model_id, flux, path, 'one flux of 3, using flux of {} m3/d'.format(flux))

        print('done running models took {} minutes'.format((t-time.time())/60))

    else:
        base_paths = []
        for flux in fluxes:
            path = os.path.join(base_dir, 'flux_{}'.format(flux))
            base_paths.append(path)

    #### extract the data ####
    outpaths = []
    for path in base_paths:
        out_path = os.path.join(data_out_dir,'sd_grid_data_{}.csv'.format(os.path.basename(path)))
        outpaths.append(out_path)
        if amalg_results:
            calc_str_dep_all_wells_grid(out_path, path)

    print('finished extracting data')
    # add metadata
    if amalg_results:
        for path in base_paths:
            save_sd_metadata(os.path.join(data_out_dir,'sd_grid_metadata_{}.csv'.format(os.path.basename(path))),path)

    # krig and plot all data
    nc_paths = []
    for path in outpaths:
        rd_path = os.path.join(os.path.dirname(path),'{}_{}'.format(model_id,os.path.basename(path)))
        nc_path = krig_plot_sd_grid(rd_path,data_out_dir)
        nc_paths.append(nc_path)

    # make shapefiles of zones
    for ncpath in nc_paths:
        outpath = os.path.join(os.path.dirname(nc_path),
                               'interpolated_zones_flux_{}.shp'.format(os.path.basename(nc_path).split('_')[-1].replace('.nc','')))
        layers = range(10)
        cutoffs = [20, 40, 60, 90]
        make_sd_zone(nc_path,layers,outpath,cutoffs)


