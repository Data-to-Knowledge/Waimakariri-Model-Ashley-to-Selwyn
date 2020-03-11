# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/10/2017 3:17 PM
"""

from __future__ import division
from waimak_extended_boundry.model_run_tools import modflow_converged
from glob import glob
import pandas as pd
import os

def save_sd_metadata(outpath, data_dir):
    """
    save the sd metatdata (E.g. convergence)
    :param outpath:
    :param data_dir:
    :return:
    """
    all_paths = glob('{}/*/*.nam'.format(data_dir))
    all_paths = [e.replace('.nam', '') for e in all_paths]
    names = [os.path.basename(e) for e in all_paths]
    model_id = os.path.basename(all_paths[0]).split('_')[0]
    status = []
    for path in all_paths:
        if os.path.getsize('{}.cbc'.format(path)) == 0:  # the model fell over and nothing was written to the cbc
            t = 'did not complete run'
        else:
            t = 'converged'
            temp = modflow_converged('{}.list'.format(path))
            if not temp:
                t = 'did not converge'
        status.append(t)

    outdata = pd.DataFrame(index=names,data=status,columns=['run_status'])
    outdata.to_csv(outpath)

