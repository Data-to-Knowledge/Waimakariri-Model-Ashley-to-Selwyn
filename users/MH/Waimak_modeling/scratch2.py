"""
Author: matth
Date Created: 27/04/2017 8:37 AM
"""

from __future__ import division

from glob import glob
from users.MH.Waimak_modeling.models.extended_boundry.nsmc_exploration_results.combine_nsmc_results.ucn_netcdf import _get_sfr_con_map
import os

if __name__ == '__main__':
    paths = glob(r"T:\Temp\temp_gw_files\sobs_check\NsmcReal000005_pc5_80.sobs")
    for path in paths:
        print os.path.basename(path)
        _get_sfr_con_map(path)
