"""
Author: matth
Date Created: 29/05/2017 8:08 AM
"""

from __future__ import division
from plot_heads import plt_transient_heads
import glob
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir
import os
plot_dir = "{}/s_explore_transient_plots".format(base_mod_dir)

glob_path = "{}/s_explore_transient/perc*\*.hds".format(base_mod_dir)

hd_paths = glob.glob(glob_path)

for hd in hd_paths:
    print('plotting for {}'.format(os.path.basename(os.path.dirname(hd))))
    plt_transient_heads(hd,'{}/{}'.format(plot_dir,os.path.basename(os.path.dirname(hd))),[0,1,2])

