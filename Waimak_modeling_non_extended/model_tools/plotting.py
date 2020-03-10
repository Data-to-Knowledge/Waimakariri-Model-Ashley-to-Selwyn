"""
Author: matth
Date Created: 10/04/2017 8:31 AM
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#from read_write import no_flow
from basic_tools import get_model_x_y
import flopy_mh as flopy
from m_wraps.base_modflow_wrapper import get_base_mf_ss


def plt_matrix(array, vmin=None, vmax=None, **kwargs):
    if vmax is None:
        vmax = np.nanmax(array)
    if vmin is None:
        vmin = np.nanmin(array)

    fig, (ax) = plt.subplots(figsize=(18.5, 9.5))
    ax.set_aspect('equal')
    model_xs, model_ys = get_model_x_y(array.shape[0],array.shape[1])
    pcm = ax.pcolormesh(model_xs, model_ys, array,
                        cmap='plasma', vmin=vmin, vmax=vmax, **kwargs)
    fig.colorbar(pcm, ax=ax, extend='max')
    ax.contour(model_xs, model_ys, no_flow)

    return fig, ax


def plt_default_map(m, lyr, a=None, plt_grid=False, vmin=None, vmax=None, title=None, masked_vals = None):
    """
    only for steady state MF at the moment?
    :param m:
    :param lyr:
    :param vmin:
    :param vmax:
    :return:
    """
    if vmax is None and a is not None:
        vmax = np.nanmax(a)
    if vmin is None and a is not None:
        vmin = np.nanmin(a)

    fig, ax = plt.subplots(1, 1, figsize=(18.5, 9.5))
    if title is not None:
        fig.suptitle(title)
    ax.set_aspect('equal')
    mmap = flopy.plot.ModelMap(model=m, layer=lyr, ax=ax, xul=1512162.53275, yul=5215083.5772)
    if a is not None:
        pcm = mmap.plot_array(a, vmin=vmin, vmax=vmax, masked_values=masked_vals, alpha=1, cmap='plasma', edgecolor='face')
        fig.colorbar(pcm, ax=ax, extend='max')
    bc_alpha = 0.25
    mmap.plot_bc('WEL', alpha=bc_alpha, color='k', edgecolors=None)
    mmap.plot_bc('STR', plotAll=True, alpha=bc_alpha, color='k', edgecolors=None)
    if plt_grid:
        mmap.plot_grid(alpha=0.25)

    mmap.plot_ibound()
    labxs = [1549031.310, 1534920.116, 1566739.476, 1572771.320, 1527836.849, 1536524.918, 1562201.758]
    labys = [5194938.736, 5206061.678, 5205784.988, 5197318.271, 5184701.203, 5183594.443, 5184977.893]
    labs = ['Cust', 'Oxf', 'Rang', 'Kai', 'Dar', 'Kir', 'Air']
    for x, y, s in zip(labxs, labys, labs):
        txt = ax.text(x, y, s, bbox=dict(facecolor='w', pad=2))
        txt.set_fontsize(8)
    fig.tight_layout()
    return fig, ax, mmap


def plt_default_xsection(m, line, a=None, c_step=25, clab=True, map_a=None, grid=False):  # set up if line is a path to shapefile # needs updating for new model
    model_xs, model_ys = get_model_x_y()
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(18.5, 9.5))
    ax.set_aspect('equal')
    mmap = flopy.plot.ModelMap(model=m, layer=0, ax=ax, xul=1512162.53275, yul=5215083.5772)
    if map_a is not None:
        pcm = mmap.plot_array(map_a, vmin=np.nanmin(map_a), vmax=np.nanmax(map_a),
                              alpha=0.5, cmap='plasma', edgecolor='face')
        fig.colorbar(pcm, ax=ax, extend='max')
    mmap.plot_bc('WEL', edgecolors=None)
    mmap.plot_bc('STR', edgecolors=None)

    mmap.plot_ibound()
    labxs = [1549031.310, 1534920.116, 1566739.476, 1572771.320, 1527836.849, 1536524.918, 1562201.758]
    labys = [5194938.736, 5206061.678, 5205784.988, 5197318.271, 5184701.203, 5183594.443, 5184977.893]
    labs = ['Cust', 'Oxf', 'Rang', 'Kai', 'Dar', 'Kir', 'Air']
    for x, y, s in zip(labxs, labys, labs):
        txt = ax.text(x, y, s, bbox=dict(facecolor='w', pad=2))
        txt.set_fontsize(8)
    if 'column' in line.keys():
        col = line.values()[0]
        ax.plot(model_xs[:,col], model_ys[:,col], color='red')

    elif 'row' in line.keys():
        row = line.values()[0]
        ax.plot(model_xs[row,:], model_ys[row,:], color='red')
    else:
        line_vals = line.values()[0]
        ax.plot(line_vals[:,0], line_vals[:,1], color='red')



    mc = flopy.plot.ModelCrossSection(ax=ax2, line=line, model=m, xul=1512162.53275, yul=5215083.5772)
    mc.plot_bc('WEL')
    mc.plot_bc('STR')
    if grid:
        mc.plot_grid()
    if a is not None:
        pcm = mc.plot_array(a, cmap='plasma', alpha=0.5)
        fig.colorbar(pcm, ax=ax, extend='max')
        levels = range(np.nanmin(a), np.nanmax(a), c_step)
        levels = [np.round(e) for e in levels]
        ct = mc.contour_array(a, colors='k', alpha=0.75, levels=levels)
        if clab:
            ax.clabel(ct, fontsize=9, inline=1)

    mc.plot_ibound()
    fig.tight_layout()

    return fig, ax, ax2, mc, mmap


if __name__ == '__main__':
    con_data_org = flopy.utils.UcnFile(r"C:\Users\MattH\Downloads\MT3D001.UCN")
    con_data = con_data_org.get_data(kstpkper=con_data_org.get_kstpkper()[-1])
    con_data[con_data > 10000] = np.nan
    m = get_base_mf_ss()
    fig, ax, mmap = plt_default_map(m, 0, con_data[0], plt_grid=False, vmax=10, vmin=0)

    print('done')
