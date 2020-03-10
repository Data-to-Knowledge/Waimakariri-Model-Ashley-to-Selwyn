"""
Author: matth
Date Created: 18/05/2017 9:02 AM
"""

from __future__ import division
import numpy as np
import flopy
import os
import shutil
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import base_mod_dir, sdp

if __name__ == '__main__': # this isn't working super well with lines debug
    hds = flopy.utils.HeadFile("{}/base_model_runs/base_ss_mf/base_SS.hds".format(sdp))
    head = hds.get_data(hds.get_kstpkper()[-1])
    head[np.isclose(head, -999)] = np.nan
    labxs = [1549031.310, 1534920.116, 1566739.476, 1572771.320, 1527836.849, 1536524.918, 1562201.758]
    labys = [5194938.736, 5206061.678, 5205784.988, 5197318.271, 5184701.203, 5183594.443, 5184977.893]
    test = np.array(zip(labxs, labys))[0:3]
    m = mt.wraps.mflow.get_base_mf_ss()
    pout_col = '{}/column_cross_sections'.format(base_mod_dir)
    if os.path.exists(pout_col):
        shutil.rmtree(pout_col)
    os.makedirs(pout_col)
    pout_row = '{}/row_cross_sections'.format(base_mod_dir)
    if os.path.exists(pout_row):
        shutil.rmtree(pout_row)
    os.makedirs(pout_row)

    for i in range(0,190,10):
        try:
            fig,ax,ax2,mc, mmap = mt.plt_default_xsection(m, {'row': i}, head, c_step=5)
            fig.savefig('{}/cross_row_{:03d}'.format(pout_row, i))
            mt.plt.close(fig)
        except:
            pass

    for i in range(10,365,10):
        try:
            fig,ax,ax2,mc, mmap = mt.plt_default_xsection(m, {'column': i}, head, c_step=5)
            fig.savefig('{}/cross_col_{:03d}'.format(pout_col, i))
            mt.plt.close(fig)
        except:
            pass

