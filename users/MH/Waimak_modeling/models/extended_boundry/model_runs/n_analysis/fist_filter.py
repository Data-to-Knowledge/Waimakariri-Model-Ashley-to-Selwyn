# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/04/2018 9:11 AM
"""

from __future__ import division
from core import env
import pandas as pd
import os
from glob import glob
import shutil

if __name__ == '__main__':
    # this didn't really do much! (too low)
    base_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_vs_wai_regressions"

    for rt in ['stream', 'well']:
        outdir = os.path.join(base_dir,"{r}s".format(r=rt),'possible_adjustments')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        path = os.path.join(base_dir, "{r}s\{r}_regression_data.csv".format(r=rt))
        data = pd.read_csv(path, index_col=0)
        data = data.loc[data.adj_r2 > 0.1]
        data.to_csv(os.path.join(outdir, "{r}_reduced_regression_data.csv".format(r=rt)))
        for site in data.index:
            base_plt = os.path.join(base_dir, "{}s/plots/{}.png".format(rt, site))
            shutil.copy(base_plt, os.path.join(outdir,'{}.png'.format(site)))
