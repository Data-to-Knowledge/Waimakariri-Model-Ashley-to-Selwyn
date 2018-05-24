# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 24/05/2018 2:22 PM
"""

from __future__ import division
from core import env
import bokeh
import numpy as np

from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Select
from bokeh.layouts import column, widgetbox, layout

from osgeo.gdal import Open
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index


def get_base_data():
    raise NotImplementedError


def plot_reduction():
    ds = Open(
        env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\topo250small.tif"))
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    maxx = gt[0] + width * gt[1] + height * gt[2]
    maxy = gt[3]

    image = np.flipud(ds.ReadAsArray())
    if image.ndim == 3:  # if a rgb image then plot as greyscale
        image = image.mean(axis=0)
    ll = (minx, miny)
    ur = (maxx, maxy)

    x_range = (1516000, 1576500)
    y_range = [5186000, 5214500]
    p = figure(tools=["pan,wheel_zoom,reset,save"], active_scroll='wheel_zoom',
               match_aspect=True)

    # must give a vector of image data for image parameter
    # plot basemap
    p.image(image=[image], x=minx, y=miny, dw=maxx - minx, dh=maxy - miny, palette="Greys256")

    y = smt.uly - smt.grid_space * smt.rows
    test = np.flipud(get_zone_array_index('waimak'))
    img = p.image(image=[test], x=smt.ulx, y=y, dw=smt.cols * smt.grid_space, dh=smt.rows * smt.grid_space,
                  palette="Plasma256",
                  alpha=0.5)
    image_ds = img.data_source

    targets = {
        'k_harpers': [6.9],
        'k_islands': [6.9],
        'ohoka': [6.9],
        'court': [6.9],
        'cust': [6.9],
        'cam': [2.4],
        'public': [5.6],
        'private': [5.6],
        'waimak': [0]}
    target_source = ColumnDataSource(targets)
    prediction_source = ColumnDataSource()  # todo think about this cannot handle any deeper objects/abstraction
    # I will need to have a use prediction souce which has keys (ohoka ect) and values (2d numpy array)
    # I will also need to have a base source for each group (waimak, ect) which has keys: 25%: array or value.
    # there will need to be a different call back function for each target select that pulls the apropriate map
    # and updates the
    # use value data source.

    def callback_tar(im_source=image_ds, tar_source=target_source, pred_source=prediction_source):  # todo
        # update the target data source with the new target
        new_tar = cb_obj.value
        name = cb_obj.name
        if name == 'waimak':
            if new_tar == 'current':
                val = 0
            else:
                val = 28
        elif name == 'public' or name == 'private':
            if new_tar == 'no target':
                val = 999
            else:
                val = float(new_tar)
        else:
            val = float(new_tar)
        tar_source.data[name] = val

        # now take the targets generate the maps and change the im source #todo make a function of this
        raise NotImplementedError

    sel_tar_kai_harp = Select(name='k_harpers', title='Kaiapoi Harpers Target', value='6.9',
                              options=['6.9', '3.8', '1.0'])
    sel_tar_kai_isl = Select(name='k_islands', title='Kaiapoi Island Target', value='6.9',
                             options=['6.9', '5.4', '3.8', '1.0'])
    sel_tar_ohoka = Select(name='ohoka', title='Ohoka Target', value='6.9', options=['6.9', '4.5', '3.8', '1.0'])
    sel_tar_cour = Select(name='court', title='Courtenay Target', value='6.9', options=['6.9', '3.8', '3.1', '1.0'])
    sel_tar_cust = Select(name='cust', title='Cust Main Drain Target', value='6.9',
                          options=['6.9', '4.7', '3.8', '1.0'])
    sel_tar_cam = Select(name='cam', title='Cam River Target', value='2.4', options=['2.4', '1.5', '1.0'])
    sel_tar_public = Select(name='public', title='Public Wells Target', value='5.65',
                            options=['no target', '10', '8.5', '5.65'])
    sel_tar_private = Select(name='private', title='Private Wells Target', value='5.65',
                             options=['no target', '8.5', '7.1', '5.65'])
    sel_tar_wai = Select(name='waimak', title='Waimakariri Target', value='current',
                         options=['current', 'periphyton control'])

    sel_likely_kai_harp = Select(name='k_harpers', title='Kaiapoi Harpers Likelyhood', value='50%',
                                 options=['5%', '50%', '95%'])
    sel_likely_kai_isl = Select(name='k_islands', title='Kaiapoi Island Likelyhood', value='50%',
                                options=['5%', '50%', '95%'])
    sel_likely_ohoka = Select(name='ohoka', title='Ohoka Likelyhood', value='50%',
                              options=['5%', '50%', '95%'])
    sel_likely_cour = Select(name='court', title='Courtenay Likelyhood', value='50%',
                             options=['5%', '50%', '95%'])
    sel_likely_cust = Select(name='cust', title='Cust Main Drain Likelyhood', value='50%',
                             options=['5%', '50%', '95%'])
    sel_likely_cam = Select(name='cam', title='Cam River Likelyhood', value='50%',
                            options=['5%', '50%', '95%'])
    sel_likely_public = Select(name='public', title='Public Wells Likelyhood', value='50%',
                               options=['5%', '50%', '95%'])
    sel_likely_private = Select(name='private', title='Private Wells Likelyhood', value='50%',
                                options=['5%', '50%', '95%'])
    layout2 = layout([p,
                      [sel_tar_kai_harp, sel_likely_kai_harp],
                      [sel_tar_kai_isl, sel_likely_kai_isl],
                      [sel_tar_ohoka, sel_likely_ohoka],
                      [sel_tar_cour, sel_likely_cour],
                      [sel_tar_cust, sel_likely_cust],
                      [sel_tar_cam, sel_likely_cam],
                      [sel_tar_public, sel_likely_public],
                      [sel_tar_private, sel_likely_private],
                      sel_tar_wai
                      ])

    # save
    #todo need legend
    output_file(r"C:\Users\MattH\Downloads\test_bokeh.html", title="image.py example")

    show(layout2)  # open a browser


if __name__ == '__main__':
    test = ColumnDataSource(
        {
            'kaiapoi': np.array([[1, 2], [3, 4]]),
            'ohoka': np.array([[1, 2], [3, 4]])
        })
    plot_reduction()
