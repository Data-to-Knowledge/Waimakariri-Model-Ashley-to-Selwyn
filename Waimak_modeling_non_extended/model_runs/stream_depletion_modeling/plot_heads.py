"""
Author: matth
Date Created: 25/05/2017 9:05 AM
"""

from __future__ import division
import users.MH.Waimak_modeling.model_tools as mt
import flopy
import numpy as np
import matplotlib as mpl
import os


def plt_transient_heads(infile, outdir, layers, stress_to_month):
    """
    plt all heads from layer one to the outdir. 
    :param infile: .hds file 
    :param outdir:  dir to save plots in
    :param layers: layers to plot either int or iterable of int
    :return: 
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    m = mt.get_base_mf_ss()
    all_heads = flopy.utils.HeadFile(infile)
    cmap = mpl.colors.ListedColormap(['1', 'w'])
    layers = np.atleast_1d(layers)
    for lay in layers:
        hds0 = all_heads.get_alldata(lay)
        dry = np.zeros(hds0.shape)
        dry_idx = np.isclose(hds0,-888.)
        hds0[dry_idx] = np.nan
        dry[dry_idx] = 1
        vmax = np.nanmax(hds0)
        for i, stp in enumerate(all_heads.get_kstpkper()):
            title = 'layer {} {} step {}'.format(lay,stress_to_month[stp[1]], stp[0])
            fig, ax, mmap = mt.plt_default_map(m,lay,hds0[i], vmax=vmax, vmin=0, title=title, masked_vals=[-999,-888])
            mmap.plot_array(dry[i], masked_values=[0], cmap=cmap)
            fig.savefig("{}/lay{:02d}_per_{:02d}_stp_{:02d}.png".format(outdir,lay,stp[1],stp[0]))
            mt.plt.close(fig)

if __name__ == '__main__':
    stress_to_month_2 = {0: 10,
                       1: 11,
                       2: 12,
                       3: 1,
                       4: 2}

    plt_transient_heads(r"C:\Users\MattH\Desktop\Waimak_modeling\python_models\well_by_well_str_dep\remove_test\remove_test.hds",
                        r"C:\Users\MattH\Desktop\Waimak_modeling\python_models\well_by_well_str_dep\remove_test_plots",
                        [0,1,2,3],
                        stress_to_month=stress_to_month_2)