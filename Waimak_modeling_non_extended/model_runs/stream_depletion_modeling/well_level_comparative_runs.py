"""
Author: matth
Date Created: 9/06/2017 10:47 AM
"""

from __future__ import division
import glob
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import results_dir, base_mod_dir
import os
import numpy as np
import warnings

def main():
    model_paths = glob.glob("{}/sd_comparative_runs/*/*.hds".format(base_mod_dir))

    wells = None # define as list
    all_wells = mt.get_all_well_data()
    well_data = all_wells[np.in1d(all_wells.index.values, wells)]
    if len(well_data.index) != len(wells):
        warnings.warn('the following wells were not found in all well data {}'.format(set(wells)-set(well_data.index.values)))
    well_hydrograph = {}
    for path in model_paths:
        well_hydrograph[os.path.basename(path)] = mt.calc_well_hydrographs(path, well_data)

    outdir = '{}/stream_depletion/stream_depletion_comparative_well_data'.format(results_dir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    stress_to_month = {0: 'steady_state',
                       1: 'jul',
                       2: 'aug',
                       3: 'sep',
                       4: 'oct',
                       5: 'nov',
                       6: 'dec',
                       7: 'jan',
                       8: 'feb',
                       9: 'mar',
                       10: 'apr',
                       11: 'may',
                       12: 'jun'}

    for key in well_hydrograph.keys():
        data = well_hydrograph[key]
        data['month'] = None
        for i in data.index:
            data.loc[i].loc['month'] = stress_to_month[i[1]]
        data.to_csv('{}/{}_all_str_flow.csv'.format(outdir,key))

if __name__ == '__main__':
    main()
