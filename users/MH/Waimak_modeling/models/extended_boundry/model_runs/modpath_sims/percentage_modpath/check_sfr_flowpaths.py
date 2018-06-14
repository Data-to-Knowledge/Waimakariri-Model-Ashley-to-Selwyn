"""
Author: matth
Date Created: 12/12/2017 3:27 PM
"""

from __future__ import division

import flopy_mh as flopy
import numpy as np
import pandas as pd

from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.modpath_wrapper import \
    export_paths_to_shapefile
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.percentage_modpath.run_emulator import \
    get_group_number

if __name__ == '__main__':
    bnd_type = np.loadtxt("D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try_bnd_type.txt")
    temp = bnd_type == 1 # well is 1 sfr is 2
    groups = get_group_number(index=temp, bd_type=bnd_type)
    end = pd.DataFrame(flopy.utils.EndpointFile(r"D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try.mpend").get_alldata())

    pid = end.loc[np.in1d(end.particlegroup,groups),'particleid'].values

    export_paths_to_shapefile(r"D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try.mppth",
                              r"T:\Temp\temp_gw_files\well_pathlines.shp",particle_ids=pid)

    bnd_type = np.loadtxt("D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try_bnd_type.txt")
    temp = bnd_type == 2 # well is 1 sfr is 2
    groups = get_group_number(index=temp, bd_type=bnd_type)
    end = pd.DataFrame(flopy.utils.EndpointFile(r"D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try.mpend").get_alldata())

    pid = end.loc[np.in1d(end.particlegroup,groups),'particleid'].values

    export_paths_to_shapefile(r"D:\mh_waimak_models\modpath_emulator\NsmcBase_first_try.mppth",
                              r"T:\Temp\temp_gw_files\sfr_pathlines.shp",particle_ids=pid)