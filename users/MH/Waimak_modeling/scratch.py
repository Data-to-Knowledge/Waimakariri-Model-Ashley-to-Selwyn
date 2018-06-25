from __future__ import division
from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results.ucn_netcdf import _get_sfr_con_map
import shutil
import os

if __name__ == '__main__':

    paths = glob("D:\mh_waimak_models\dairy_lowdrain\*\*.sobs")
    outdir = os.path.join(r"T:\Temp\temp_gw_files", 'sobs_check')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for path in paths:
        shutil.copyfile(path, os.path.join(outdir,os.path.basename(path)))
    print('done')
