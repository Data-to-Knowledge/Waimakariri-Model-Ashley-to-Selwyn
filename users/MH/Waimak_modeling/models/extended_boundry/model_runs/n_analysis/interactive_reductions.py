# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 24/05/2018 2:22 PM
"""

from __future__ import division
from core import env
import numpy as np
from glob import glob
import os

from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models.widgets import Select
from bokeh.layouts import column, widgetbox, layout

from osgeo.gdal import Open
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
from get_current_pathway_n import get_mt3d_current_pathway_n
from percentage_reduction_maps import wdc_wells, private_wells

paths_to_shp = {
    'k_harpers': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\kaiapoi_harpers_s.shp",
    'k_islands': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\kaiapoi_island_s.shp",
    'ohoka': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\ohoka_island_s.shp",
    'court': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\courtenay_kaiapoi_s.shp",
    'cust': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\cust_skewbridge.shp",
    'cam': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\second_tranche\cam_marshes_s.shp",
    'public': glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\wdc_use_mix\*.shp"),
    'private': glob(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\private_wells_90_named_right\*.shp"),
    'waimak': r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\waimakariri_river\waimakariri.shp"
}


def get_base_data():
    areas = ['k_harpers',
             'k_islands',
             'ohoka',
             'court',
             'cust',
             'cam', ]
    pc5pas = [True, False]
    likelihoods = ['5%', '50%', '95%']
    outdata = {}
    for pc5 in pc5pas:
        pa = 'pc5' if pc5 else 'wopc5'
        current_pathway = get_mt3d_current_pathway_n(add_pc5=pc5)
        for like in likelihoods:
            outdata['k_harpers_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['k_harpers'], 'Id', True))) * current_pathway.loc[
                                                              'kaiapoi_harpers_s', like]
            outdata['k_islands_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['k_islands'], 'Id', True))) * current_pathway.loc[
                                                              'kaiapoi_island_s', like]
            outdata['ohoka_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['ohoka'], 'Id', True))) * current_pathway.loc[
                                                          'ohoka_island_s', like]
            outdata['court_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['court'], 'Id', True))) * current_pathway.loc[
                                                          'courtenay_kaiapoi_s', like]
            outdata['cust_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['cust'], 'Id', True))) * current_pathway.loc[
                                                         'cust_skewbridge', like]
            outdata['cam_{}_{}'.format(like, pa)] = np.flipud(np.isfinite(
                smt.shape_file_to_model_array(paths_to_shp['cam'], 'Id', True))) * current_pathway.loc[
                                                        'cam_bramleys_s', like]

            # add private_wells
            private_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\private_wells_90_named_right"
            private_reductions = []
            for well in private_wells:
                temp = smt.shape_file_to_model_array(os.path.join(private_dir, well + '.shp'),
                                                     'Id', True)
                temp[np.isfinite(temp)] = current_pathway.loc[well, like]
                private_reductions.append(temp[np.newaxis])
            private_reductions = np.flipud(np.nanmax(np.concatenate(private_reductions, axis=0), axis=0))
            outdata['private_{}_{}'.format(like, pa)] = private_reductions

            # add public_wells
            wdc_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\wdc_use_mix"
            wdc_reductions = []
            for well in wdc_wells:
                temp = smt.shape_file_to_model_array(os.path.join(wdc_dir, well.replace('wdc_', '') + '.shp'),
                                                     'Id', True)
                temp[np.isfinite(temp)] = current_pathway.loc[well, like]
                wdc_reductions.append(temp[np.newaxis])
            wdc_reductions = np.flipud(np.nanmax(np.concatenate(wdc_reductions, axis=0), axis=0))
            outdata['public_{}_{}'.format(like, pa)] = wdc_reductions

    outdata['waimak'] = np.flipud(np.isfinite(smt.shape_file_to_model_array(paths_to_shp['waimak'], 'Id', True)).astype(float))
    return outdata


def create_reduction_reason(tar_ds, pred_ds, version='reduction'):
    areas = ['k_harpers',
             'k_islands',
             'ohoka',
             'court',
             'cust',
             'cam',
             'public',
             'private']

    outdata = []
    for area in areas:
        likelihood = tar_ds.data['lik_{}'.format(area)][0]
        pa = 'pc5' if tar_ds.data['pc5pa'.format(area)][0] else 'wopc5'
        target = tar_ds.data['tar_{}'.format(area)][0]
        org = pred_ds.data['{}_{}_{}'.format(area, likelihood, pa)]
        org = (1 - target / org) * 100
        org[org < 0] = 0
        outdata.append(org[np.newaxis])
    org = pred_ds.data['waimak'].astype(float)
    org[np.isclose(org, 0)] = 0
    outdata.append((org * tar_ds.data['tar_waimak'][0])[np.newaxis])
    if version == 'reduction':
        return np.nanmax(np.concatenate(outdata, axis=0), axis=0)
    else:
        raise NotImplementedError('version: {}'.format(version))


if __name__ == '__main__':
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
    # todo add hover?
    # must give a vector of image data for image parameter
    # plot basemap
    p.image(image=[image], x=minx, y=miny, dw=maxx - minx, dh=maxy - miny, palette="Greys256")

    targets = {
        # targets
        'tar_k_harpers': [6.9],
        'tar_k_islands': [6.9],
        'tar_ohoka': [6.9],
        'tar_court': [6.9],
        'tar_cust': [6.9],
        'tar_cam': [2.4],
        'tar_public': [999],
        'tar_private': [999],
        'tar_waimak': [0],
        # pc5pa
        'pc5pa': [True],

        # liklyhoods
        'lik_k_harpers': ['50%'],
        'lik_k_islands': ['50%'],
        'lik_ohoka': ['50%'],
        'lik_court': ['50%'],
        'lik_cust': ['50%'],
        'lik_cam': ['50%'],
        'lik_public': ['50%'],
        'lik_private': ['50%']

    }
    target_source = ColumnDataSource(targets)
    prediction_source = ColumnDataSource(
        get_base_data())  # todo think about this cannot handle any deeper objects/abstraction

    # I will need to have a use prediction souce which has keys (ohoka ect) and values (2d numpy array)
    # I will also need to have a base source for each group (waimak, ect) which has keys: 25%: array or value.
    # there will need to be a different call back function for each target select that pulls the apropriate map
    # and updates the
    # use value data source.
    # other thoughts simply make a static dictionary of prediction arrays under pc5pa and likelhoods
    # keep a register of boolean pa and which percentage level for each receptor in the target_source dictionary


    sel_pa = Select(name='pc5pa', title='PC5PA rules', value='With PC5PA',
                    options=['With PC5PA', 'Without PC5PA'])
    sel_tar_kai_harp = Select(name='tar_k_harpers', title='Kaiapoi Harpers Target', value='6.9',
                              options=['6.9', '3.8', '1.0'])
    sel_tar_kai_isl = Select(name='tar_k_islands', title='Kaiapoi Island Target', value='6.9',
                             options=['6.9', '5.4', '3.8', '1.0'])
    sel_tar_ohoka = Select(name='tar_ohoka', title='Ohoka Target', value='6.9', options=['6.9', '4.5', '3.8', '1.0'])
    sel_tar_cour = Select(name='tar_court', title='Courtenay Target', value='6.9', options=['6.9', '3.8', '3.1', '1.0'])
    sel_tar_cust = Select(name='tar_cust', title='Cust Main Drain Target', value='6.9',
                          options=['6.9', '4.7', '3.8', '1.0'])
    sel_tar_cam = Select(name='tar_cam', title='Cam River Target', value='2.4', options=['2.4', '1.5', '1.0'])
    sel_tar_public = Select(name='tar_public', title='Public Wells Target', value='no target',
                            options=['no target', '10', '8.5', '5.65'])
    sel_tar_private = Select(name='tar_private', title='Private Wells Target', value='no target',
                             options=['no target', '8.5', '7.1', '5.65'])
    sel_tar_wai = Select(name='tar_waimak', title='Waimakariri Target', value='current',
                         options=['current', 'periphyton control'])

    sel_likely_kai_harp = Select(name='lik_k_harpers', title='Kaiapoi Harpers Likelyhood', value='50%',
                                 options=['5%', '50%', '95%'])
    sel_likely_kai_isl = Select(name='lik_k_islands', title='Kaiapoi Island Likelyhood', value='50%',
                                options=['5%', '50%', '95%'])
    sel_likely_ohoka = Select(name='lik_ohoka', title='Ohoka Likelyhood', value='50%',
                              options=['5%', '50%', '95%'])
    sel_likely_cour = Select(name='lik_court', title='Courtenay Likelyhood', value='50%',
                             options=['5%', '50%', '95%'])
    sel_likely_cust = Select(name='lik_cust', title='Cust Main Drain Likelyhood', value='50%',
                             options=['5%', '50%', '95%'])
    sel_likely_cam = Select(name='lik_cam', title='Cam River Likelyhood', value='50%',
                            options=['5%', '50%', '95%'])
    sel_likely_public = Select(name='lik_public', title='Public Wells Likelyhood', value='50%',
                               options=['5%', '50%', '95%'])
    sel_likely_private = Select(name='lik_private', title='Private Wells Likelyhood', value='50%',
                                options=['5%', '50%', '95%'])

    # plot reduction
    y = smt.uly - smt.grid_space * smt.rows
    test = create_reduction_reason(target_source, prediction_source, version='reduction')
    img = p.image(image=[test], x=smt.ulx, y=y, dw=smt.cols * smt.grid_space, dh=smt.rows * smt.grid_space,
                  palette="Plasma256",
                  alpha=0.5)  # todo try to figure out different plotting levels, and how to set a datasource for reason
    image_ds = img.data_source


    def callback_tar(source=image_ds, window=None):  # todo I may need to pass the
        # update the target data source with the new target
        #new_tar = cb_obj.value
        #name = cb_obj.name
        #if 'lik' in name:  # pass the likelihoods straight across
        #    val = new_tar
        #elif name == 'tar_waimak':
        #    if new_tar == 'current':
        #        val = 0
        #    else:
        #        val = 28
        #elif name == 'tar_public' or name == 'tar_private':
        #    if new_tar == 'no target':
        #        val = 999
        #    else:
        #        val = float(new_tar)
        #elif name == 'pc5pa':
        #    if new_tar == 'With PC5PA':
        #        val = True
        #    else:
        #        val = False
        #else:
        #    val = float(new_tar)
        #tar_source.data[name] = val
#
        #areas = ['k_harpers',
        #         'k_islands',
        #         'ohoka',
        #         'court',
        #         'cust',
        #         'cam',
        #         'public',
        #         'private']
#
        #outdata = []
        #for area in areas:
        #    likelihood = tar_source.data['lik_{}'.format(area)][0]
        #    pa = 'pc5' if tar_source.data['pc5pa'.format(area)][0] else 'wopc5'
        #    target = tar_source.data['tar_{}'.format(area)][0]
        #    org = pred_source.data['{}_{}_{}'.format(area, likelihood, pa)]
        #    for i in range(364):
        #        for j in range(365):
        #            org[i,j] = (1 - target / org[i,j]) * 100
        #            if org[i,j] < 0:
        #                org[i,j] = 0
        #    outdata.append(org)
        #org = pred_source.data['waimak'] #todo make sure this is not a list
        #for i in range(364):
        #    for j in range(365):
        #        org[i, j] = org[i,j] * tar_source.data['tar_waimak'][0]
        #outdata.append(org)
        # temp = im_source.data['image'][0]
        temp2 = []
        for i in range(364):
            temp3 = []
            for j in range(365):
                temp3.append(i)
            temp2.append(temp3)
        source.data.image = temp2 # this breaks it for some reason, but you see a change
                #source.data['image'][i][j] = 0 #try new not mutating
                #max(outdata[0][i][j],outdata[1][i][j],outdata[2][i][j],outdata[3][i][j],outdata[4][i][j],
                #outdata[5][i][j],outdata[6][i][j]) #todo check how many arrays
        source.change.emit()
    #todo put everything in a singel data source
    # todo add the update on change thing
    sel_tar_private.js_on_change('value', CustomJS.from_py_func(callback_tar))

    layout2 = layout([p,
                      sel_pa,
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
    # todo need legend
    output_file(r"C:\Users\MattH\Downloads\test_bokeh.html", title="image.py example")

    show(layout2)  # open a browser
