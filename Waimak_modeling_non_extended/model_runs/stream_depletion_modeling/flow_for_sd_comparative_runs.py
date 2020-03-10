"""
Author: matth
Date Created: 8/06/2017 3:43 PM
"""

from __future__ import division
import glob
import users.MH.Waimak_modeling.model_tools as mt
import pandas as pd
from users.MH.Waimak_modeling.supporting_data_path import results_dir, base_mod_dir
import os
import numpy as np

def main():
    model_paths = glob.glob("{}/sd_comparative_runs/*/*.hds".format(base_mod_dir))
    model_paths = [e.replace('.hds','') for e in model_paths]

    str_points = mt.get_str_sample_points()
    idx = ['hyg' in e for e in str_points]
    str_points = list(np.array(str_points)[np.array(idx)])

    drn_points = np.array(mt.get_drn_samp_pts())
    idx = np.array(['swaz' not in e for e in drn_points])
    drn_points = list(drn_points[idx])

    stream_flows = {}
    for path in model_paths:
        stream_flows[os.path.basename(path)] = mt.streamflow_for_kskps(path,drn_points=drn_points,str_points=str_points)

    outdir = '{}/stream_depletion/stream_depletion_comparative'.format(results_dir)
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

    outdata = pd.DataFrame()
    outdata2 = pd.DataFrame()
    for key in stream_flows.keys():
        data = stream_flows[key]
        data *= 0.0115741 # convert to l/s
        data['month'] = None
        for i in data.index:
            data.loc[i].loc['month'] = stress_to_month[i[1]]
        data.to_csv('{}/{}_all_str_flow.csv'.format(outdir,key))

        temp = data.loc[1,9] #grab data at the end of march
        temp2 = data.loc[0,0]
        outdata.loc[:,key] = temp
        outdata2.loc[:,key] = temp2
    outdata = outdata.drop('month', axis=0)
    outdata2 = outdata2.drop('month', axis=0)
    outdata['gmp_percent_cmp'] = outdata['m_strong_vert-gmp_simulation']/outdata['m_strong_vert-current_abstraction']*100
    outdata2['gmp_percent_cmp'] = outdata2['m_strong_vert-gmp_simulation']/outdata2['m_strong_vert-current_abstraction']*100
    outdata['nat_percent_cmp'] = outdata['m_strong_vert-fully_nat_without_races']/outdata['m_strong_vert-current_abstraction']*100
    outdata2['nat_percent_cmp'] = outdata2['m_strong_vert-fully_nat_without_races']/outdata2['m_strong_vert-current_abstraction']*100
    outdata.index.name = 'l/s'
    outdata2.index.name = 'l/s'
    outdata = outdata.transpose()
    outdata2 = outdata2.transpose()
    outdata.to_csv('{}/all_flows_end_march.csv'.format(outdir))
    outdata2.to_csv('{}/all_flows_stead_state.csv'.format(outdir))

if __name__ == '__main__':
    main()
