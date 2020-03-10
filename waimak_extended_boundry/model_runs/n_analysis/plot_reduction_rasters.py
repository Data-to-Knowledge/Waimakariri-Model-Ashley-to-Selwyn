# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 16/05/2018 9:13 AM
"""

from __future__ import division
from osgeo.gdal import Open
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from waimak_extended_boundry import smt

def plot_reason_raster(path, outpath=None, title=None, include_interzone=True):
    data = Open(path).ReadAsArray()

    fig, (ax) = plt.subplots(figsize=(18.5, 9.5))
    if title is not None:
        ax.set_title(title)
    ax.set_aspect('equal')

    model_xs, model_ys = smt.get_model_x_y()
    ds = Open(smt.base_map_path)
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    maxx = gt[0] + width * gt[1] + height * gt[2]
    maxy = gt[3]

    image = ds.ReadAsArray()
    if image.ndim == 3:  # if a rgb image then plot as greyscale
        image = image.mean(axis=0)
    ll = (minx, miny)
    ur = (maxx, maxy)

    ax.imshow(image, extent=[ll[0], ur[0], ll[1], ur[1]], cmap='gray')  #
    if include_interzone:
        levels = np.arange(0.5,7.5,1)
    else:
        levels = np.arange(1.5,7.5,1)
    cf = ax.contourf(model_xs, model_ys, data, alpha=0.5, levels=levels, cmap='tab10')
    ticks = (np.array(levels)+ 0.5)[:-1]
    cbar = fig.colorbar(cf, ticks=ticks)
    if include_interzone:
        labs = ['Interzone', 'Wdc wells', 'Private wells', 'Springfed streams', 'Waimakariri R.', 'All equal']
    else:
        labs = ['Wdc wells', 'Private wells', 'Springfed streams', 'Waimakariri R.', 'All equal']

    cbar.ax.set_yticklabels(labs)
    cbar.ax.tick_params(labelsize=20)
    xlim,ylim = smt._get_xlim_ylim()
    ax.set_ylim([5186000, 5214500])
    ax.set_xlim((1516000,1576500))
    fig.tight_layout()

    if outpath is not None:
        fig.savefig(outpath)
    return fig,ax


def plot_reduction_raster(path, outpath=None, title=None):
    data = Open(path).ReadAsArray()

    fig, (ax) = plt.subplots(figsize=(18.5, 9.5))
    if title is not None:
        ax.set_title(title)
    ax.set_aspect('equal')

    model_xs, model_ys = smt.get_model_x_y()
    ds = Open(smt.base_map_path)
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    maxx = gt[0] + width * gt[1] + height * gt[2]
    maxy = gt[3]

    image = ds.ReadAsArray()
    if image.ndim == 3:  # if a rgb image then plot as greyscale
        image = image.mean(axis=0)
    ll = (minx, miny)
    ur = (maxx, maxy)

    ax.imshow(image, extent=[ll[0], ur[0], ll[1], ur[1]], cmap='gray')  #
    levels = range(1, 70, 10)
    cf = ax.contourf(model_xs, model_ys, data, alpha=0.5, levels=levels, cmap='viridis',extend='max')
    ticks = np.array(levels)+5
    ticks[-1] = levels[-1]+1
    cbar = fig.colorbar(cf, ticks=ticks)
    labs = ['{}%-{}%'.format(e,e+10) for e in levels]
    labs[-1] = '>{}%'.format(levels[-1])
    cbar.ax.set_yticklabels(labs)
    cbar.ax.tick_params(labelsize=20)
    xlim,ylim = smt._get_xlim_ylim()
    ax.set_ylim([5186000, 5214500])
    ax.set_xlim((1516000,1576500))
    fig.tight_layout()

    if outpath is not None:
        fig.savefig(outpath)
    return fig,ax



def plot_all_rasters(basedir, outdir, interzone_in_reason=False):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    paths = glob(os.path.join(basedir, '*reduction.tif'))
    for path in paths:
        fig, ax = plot_reduction_raster(path,
                                        outpath=os.path.join(outdir, os.path.basename(path).replace('.tif', '.png')),
                                        title=None)
        plt.close(fig)
    paths = glob(os.path.join(basedir, '*reason.tif'))
    for path in paths:
        fig, ax = plot_reason_raster(path,
                                        outpath=os.path.join(outdir, os.path.basename(path).replace('.tif', '.png')),
                                        title=None,
                                     include_interzone=interzone_in_reason)
        plt.close(fig)



if __name__ == '__main__':
    plot_all_rasters(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_red_mt3d",
                     r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_red_mt3d\plots")
    print('done')
