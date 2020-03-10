"""
Author: matth
Date Created: 17/05/2017 8:06 AM
"""

from __future__ import division
import glob
import flopy
import users.MH.Waimak_modeling.model_tools as mt
import os
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
import shutil
import matplotlib.pyplot as plt
import numpy as np
import socket

def main():
    _comp = socket.gethostname()
    if _comp == 'HP1639':
        root_output_dir = "{}/mt3d_relaxed_vert".format(base_mod_dir)

    elif _comp == 'DHI-Runs02':
        root_output_dir = '{}/component_con_n_zero_start'.format(base_mod_dir)
    else:
        raise ValueError('unexpected computer {}'.format(socket.gethostname()))
    if not os.path.exists(root_output_dir):
        raise ValueError('output not present perhaps on VM')

    path_list = glob.glob("{}/*/*.UCN".format(root_output_dir))  #get all UCN files
    if _comp == 'DHI-Runs02':
        path_list.append("{}/first_mt3d2/cmp_layer/MT3D001.UCN".format(base_mod_dir))

    plot_output_dir = '{}/component_con_n_zero_start_plotting'.format(base_mod_dir)
    if os.path.exists(plot_output_dir):
        shutil.rmtree(plot_output_dir)  # remove the plotting dir so no versioning problems
    os.makedirs(plot_output_dir)

    m = mt.wraps.mflow.get_base_mf_ss()

    for path in path_list:
        run_name = os.path.basename(os.path.dirname(path))
        print('creating plots for {}'.format(run_name))
        out = '{}/{}'.format(plot_output_dir, run_name)
        os.makedirs(out)

        ucn = flopy.utils.UcnFile(path)
        con_data = ucn.get_data(ucn.get_kstpkper()[-1])  # get con data for the last time step
        con_data[np.isclose(con_data,1e30)] = np.nan
        con_data[con_data < -100] = np.nan

        for i in range(17):
            fig, ax, mmap = mt.plt_default_map(m,i,con_data,vmin=0,vmax=10,title='{} Layer {:02d}'.format(run_name, i+1))
            fig.savefig('{}/{}_layer{:02d}.png'.format(out,run_name,i+1))
            plt.close()
if __name__ == '__main__':
    main()