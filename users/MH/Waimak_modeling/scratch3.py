"""
Author: matth
Date Created: 22/05/2017 1:47 PM
"""

from __future__ import division
from core import env
import model_tools as mt
import flopy
import numpy as np
import matplotlib as mpl
from pykrige.ok import OrdinaryKriging




m = mt.get_base_mf_ss()
all_heads = flopy.utils.HeadFile(r"D:\mattH\python_wm_runs\test_transient\test_10_stp\test_10_stp.hds")
hds0 = all_heads.get_alldata(0)
dry = np.zeros(hds0.shape)
dry_idx = np.isclose(hds0,-888.)
hds0[dry_idx] = np.nan
dry[dry_idx] = 1
vmin = np.nanmin(hds0)
vmax = np.nanmax(hds0)
stress_to_month = {0: 'ss',
                       1: 'Jul',
                       2: 'Aug',
                       3: 'Sep',
                       4: 'Oct',
                       5: 'Nov',
                       6: 'Dec',
                       7: 'Jan',
                       8: 'Feb',
                       9: 'Mar',
                       10: 'Apr',
                       11: 'May',
                       12: 'Jun'}
cmap = mpl.colors.ListedColormap(['1', 'w'])
for i, stp in enumerate(all_heads.get_kstpkper()):
    title = '{} step {}'.format(stress_to_month[stp[1]], stp[0])
    fig, ax, mmap = mt.plt_default_map(m,0,hds0[i], vmax=vmax, vmin=0, title=title, masked_vals=[-999,-888])
    mmap.plot_array(dry[i], masked_values=[0], cmap=cmap)
    fig.savefig(r"D:\mattH\python_wm_runs\test_transient\step_10_plots\per_{}_stp_{}.png".format(stp[1],stp[0]))
    mt.plt.close(fig)
print('done')



