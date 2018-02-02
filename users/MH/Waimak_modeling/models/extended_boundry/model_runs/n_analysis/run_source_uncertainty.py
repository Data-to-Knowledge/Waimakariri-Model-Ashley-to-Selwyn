# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 31/01/2018 2:02 PM

headers for kate's data:  dairy_l,      # Dairy_light
	                      dairy_m,      # Dairy_medium
	                      dairy_pd,     # Dairy_pd
	                      sbd_l,        # Sheep_Beef_Deer_exhill_light
	                      sbd_F,        # Sheep_Beef_Deer_exhill_F
	                      sbd_S,        # Sheep_Beef_Deer_Hill_S
	                      lifestyle,    # Lifestyle
	                      doc           # Forest_DOC 
above are also the unique identifiers for the shapefile's classes with name: n_class


"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import geopandas as gpd
from glob import glob
import time
import os


def _make_dummy_file(outpath):
    headers = ['dairy_l',
               'dairy_m',
               'dairy_pd',
               'sbd_l',
               'sbd_F',
               'sbd_S',
               'lifestyle',
               'doc']
    data = {}
    for key in headers:
        data[key] = np.random.uniform(low=-100, high=200, size=10000)
    data = pd.DataFrame(data)
    data.to_csv(outpath)


def _make_dummy_shape_file():
    headers = ['dairy_l',
               'dairy_m',
               'dairy_pd',
               'sbd_l',
               'sbd_F',
               'sbd_S',
               'lifestyle',
               'doc']
    org_shape = gpd.read_file(r"C:\Users\MattH\Downloads\dummy_nload.shp")
    num = len(org_shape)
    org_shape.loc[:, 'n_class'] = [headers[e] for e in np.random.randint(0, len(headers), num)]
    org_shape.to_file(r"C:\Users\MattH\Downloads\dummy_nload_w_nclass.shp")


def calc_n_for_zone(n_zone_shp, source_area_shp_path, sims, n_load_name, outpath):
    """

    :param n_zone_shp: a geopandas object for the N zones
    :param source_area_shp_path: path to the area shapefile
    :param sims: the pre-processed data from kate
    :param outpath: path to save the raw n distribution
    :return:
    """
    if not os.path.exists(os.path.dirname(outpath)):
        os.makedirs(os.path.dirname(outpath))
    # todo write some serious assertion errors?
    zone = gpd.read_file(source_area_shp_path)
    assert len(zone) == 1, 'must have only one zone per shapefile'
    org_n_load_name = 'nload_cmp'

    # do intersection
    print('starting intersection')
    t = time.time()
    geometry = zone['geometry'].iloc[0]
    sindex = n_zone_shp.sindex
    possible_matches_index = list(sindex.intersection(geometry.bounds))
    possible_matches = n_zone_shp.iloc[possible_matches_index]
    print('took {} s to find possible matches'.format(time.time() - t))
    t = time.time()
    interest_area = gpd.overlay(possible_matches, zone,
                                how='intersection')
    # this was really slow... do a r-tree spatial: http://geoffboeing.com/2016/10/r-tree-spatial-index-python/
    print('took {} s to find specific matches'.format(time.time() - t))
    # shift to numpy
    area = interest_area.geometry.area.values
    n_type = interest_area.loc[:, 'n_class'].values
    org_nload = interest_area.loc[:, org_n_load_name] * area
    stocastic_modifiers = np.atleast_2d([sims[e] for e in (n_type)])
    temp_n_mods = np.repeat(interest_area.loc[:, n_load_name].values[:, np.newaxis],
                            stocastic_modifiers.shape[1], axis=1) * area[:,np.newaxis]
    temp_n_mods += temp_n_mods * stocastic_modifiers / 100
    n_mods = np.nansum(temp_n_mods, axis=0) / np.nansum(org_nload)
    if outpath is not None:
        np.savetxt(outpath, n_mods)
    return n_mods #todo how to check this???


def calc_all_ns(n_load_name, outdir):
    source_zone_dir = r"C:\Users\MattH\Downloads\dummy_source_zones"  # todo this is just a dummy
    n_load_path = r"C:\Users\MattH\Downloads\dummy_nload_w_nclass.shp"  # todo this is just a dummy
    sims = pd.read_csv(r"C:\Users\MattH\Downloads\dummy.csv")  # todo (from kate or dummy for now)

    sims = sims.to_dict(orient='list')
    percentiles = [0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99]
    source_paths = glob(os.path.join(source_zone_dir, '*.shp'))
    names = [os.path.basename(path).replace('.shp', '') for path in source_paths]
    n_load_layer = gpd.read_file(n_load_path)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    temp = pd.DataFrame(np.random.rand(10)).describe(percentiles=percentiles)
    outdata = pd.DataFrame(index=names, columns=temp.index)
    for path, name in zip(source_paths, names):
        print('calculating N modifiers for {}'.format(name))
        temp_n = calc_n_for_zone(n_load_layer, path, sims, n_load_name,
                                 os.path.join(outdir, 'raw_data', '{}.txt'.format(name)))
        # add the data to a summary sheet with u, sd, and some percentiles to make a quick overview of the PDF
        outdata.loc[name] = pd.Series(temp_n).describe(percentiles=percentiles)
    outdata.to_csv(os.path.join(outdir, 'n_summary_data.csv'))
    # todo what is going on with GMP? and/or other future pathways
    # todo only export modifiers
    # todo who is making the main shapefiles?

if __name__ == '__main__':
    calc_all_ns('nload_cmp', r"C:\Users\MattH\Downloads\test_n_transform")
