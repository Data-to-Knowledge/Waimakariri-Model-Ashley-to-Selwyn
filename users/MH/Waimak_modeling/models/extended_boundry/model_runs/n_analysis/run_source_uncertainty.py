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
import socket

assert socket.gethostname() == 'RDSProd03', 'must be run on RDSProd03'
import sys

repository_path = 'D:/git_repositories/matth/Ecan.Science.Python.Base'
if not repository_path in sys.path:
    sys.path.append(repository_path)
from core import env
import numpy as np
import pandas as pd
import geopandas as gpd
from glob import glob
import time
import os
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.nitrate_at_key_receptors import \
    get_well_ids, get_str_ids


def _make_dummy_file(outpath):
    """
    make a dummy file for the data Kate is going to provide
    :param outpath:
    :return:
    """
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
    """
    make a dummy shapefile for the data ogi is going to provide
    :return:
    """
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
    generate n modifiers for stocastic n
    :param n_zone_shp: a geopandas object for the N zones
    :param source_area_shp_path: path to the area shapefile
    :param sims: the pre-processed data from kate
    :param outpath: path to save the raw n distribution
    :return:
    """
    if not os.path.exists(os.path.dirname(outpath)):
        os.makedirs(os.path.dirname(outpath))
    zone = gpd.read_file(source_area_shp_path)
    assert len(zone) == 1, 'must have only one zone per shapefile'
    org_n_load_name = 'nload_cmp'

    # do intersection
    print('starting intersection')
    t = time.time()
    interest_area = spatial_overlays(n_zone_shp, zone)
    print('took {} s for intersection'.format(time.time() - t))
    # shift to numpy
    area = interest_area.geometry.area.values
    n_type = interest_area.loc[:, 'n_class'].values
    org_nload = interest_area.loc[:, org_n_load_name] * area
    stocastic_modifiers = np.atleast_2d([sims[e] for e in (n_type)])
    temp_n_mods = np.repeat(interest_area.loc[:, n_load_name].values[:, np.newaxis],
                            stocastic_modifiers.shape[1], axis=1) * area[:, np.newaxis]
    temp_n_mods += temp_n_mods * stocastic_modifiers / 100
    n_mods = np.nansum(temp_n_mods, axis=0) / np.nansum(org_nload)
    if outpath is not None:
        np.savetxt(outpath, n_mods)
    return n_mods


def calc_all_ns(sims_org, n_load_name, outdir, source_zone_dir):
    """

    :param sims_org: the simulation data (from Kate) (pd.DataFrame
    :param n_load_name: the name of the load variable to use in the shapefile
    :param outdir: the directory to save teh data
    :param source_zone_dir: the directory with all the source zone shape files.  these must be 1 shape files
    :return:
    """
    headers = ['dairy_l',
               'dairy_m',
               'dairy_pd',
               'sbd_l',
               'sbd_F',
               'sbd_S',
               'lifestyle',
               'doc']
    n_load_path = env.sci(
        'Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')

    sims = sims_org.to_dict(orient='list')
    assert set(headers) == set(sims.keys()), 'unexpected keys for sims: {} only expected: {}'.format(
        set(sims.keys()) - set(headers), headers)
    percentiles = [0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99]
    source_paths = glob(os.path.join(source_zone_dir, '*.shp'))
    names = [os.path.basename(path).replace('.shp', '') for path in source_paths]
    str_ids = get_str_ids()
    well_ids = get_well_ids()
    expected_names = np.array(
        list((set(str_ids) | set(well_ids.Zone) | set(well_ids.Zone_1) | set(well_ids.zone_2)) - {np.nan}))
    exists = np.in1d(names, expected_names)
    assert exists.all(), 'unexpected shapefile names: {} only the following allowed: {}'.format(
        np.array(names)[~exists], expected_names)
    n_load_layer = gpd.read_file(n_load_path)
    assert 'n_class' in n_load_layer.keys(), 'n_class needed in the n load layer'

    assert np.in1d(n_load_layer.n_class,
                   headers).all(), 'unexpected entries for n load layer: {} only expected {}'.format(
        set(n_load_layer.n_class) - set(headers), headers)
    assert n_load_name in n_load_layer, 'n_load_name: {} not found in n load layer, only availible: {}'.format(
        n_load_name, n_load_layer.keys())
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
    outdata.to_csv(os.path.join(outdir, 'n_modifier_summary_data.csv'))


def output_actual_n_vals(outdir, mod_dir):
    """
    load the n modifiers, and applies them to the raw n data (calculated via Nitrate at key receptors) for both the
    stocastic set and AshOpt
    :param outdir: directory to save the output
    :param mod_dir: the direcotry that holds teh n_mods (generated by calc_all_ns)
    :return:
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if os.path.exists(os.path.join(outdir, 'missing_sites.txt')):
        os.remove(os.path.join(outdir, 'missing_sites.txt'))
    # stocastic set
    # take raw data both N and modifier and put out values
    well_ids = get_well_ids()
    base_well_n = pd.read_csv(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\raw_stocastic_set_well_data.csv"),
        index_col=0).transpose()
    # well groups
    zone_sets = [set(well_ids.Zone[well_ids.Zone.notnull()]),
                 set(well_ids.Zone_1[well_ids.Zone_1.notnull()]),
                 set(well_ids.zone_2[well_ids.zone_2.notnull()])]
    outdata = {}
    for zone_set, key in zip(zone_sets, ['Zone', 'Zone_1', 'zone_2']):
        for zone in zone_set:
            try:
                modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(zone)))
            except IOError as val:
                with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                    f.write('missing: {}, exception: {}\n'.format(zone, val))
                    continue
            idxs = well_ids.loc[well_ids[key] == zone].index
            data = base_well_n.loc[idxs].mean().values
            all_n = data[:, np.newaxis] * modifiers[np.newaxis, :]
            outdata[zone] = _np_describe(all_n)
    outdata = pd.DataFrame(outdata).transpose()
    outdata.to_csv(os.path.join(outdir, 'full_stocastic_n_wells.csv'))

    str_ids = get_str_ids()
    # sfr_groups
    base_str_n = pd.read_csv(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\raw_stocastic_set_str_data.csv"),
        index_col=0)
    outdata = {}
    for str_id in str_ids:
        try:
            modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(str_id)))
        except IOError as val:
            with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                f.write('missing: {}, exception: {}\n'.format(str_id, val))
                continue
        data = base_str_n.loc[:, str_id].values
        all_n = data[:, np.newaxis] * modifiers[np.newaxis, :]
        outdata[str_id] = _np_describe(all_n)
    outdata = pd.DataFrame(outdata).transpose()
    outdata.to_csv(os.path.join(outdir, 'full_stocastic_n_strs.csv'))

    # ashopt
    # wells
    base_well_n = pd.read_csv(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\AshOpt_grouped_well_data.csv"),
        index_col=0, names=['n_con'])

    outdata = {}
    for zone_set, key in zip(zone_sets, ['Zone', 'Zone_1', 'zone_2']):
        for zone in zone_set:
            try:
                modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(zone)))
            except IOError as val:
                with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                    f.write('missing: {}, exception: {}\n'.format(zone, val))
                    continue
            data = np.atleast_1d(base_well_n.loc[zone])
            all_n = data[:, np.newaxis] * modifiers[np.newaxis, :]
            outdata[zone] = _np_describe(all_n)
    outdata = pd.DataFrame(outdata).transpose()
    outdata.to_csv(os.path.join(outdir, 'Ash_Opt_stocastic_n_wells.csv'))

    # streams
    base_str_n = pd.read_csv(env.sci(
        r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\AshOpt_stream_data.csv"),
        index_col=0)
    outdata = {}
    for str_id in str_ids:
        try:
            modifiers = np.loadtxt(os.path.join(mod_dir, 'raw_data', '{}.txt'.format(str_id)))
        except IOError as val:
            with open(os.path.join(outdir, 'missing_sites.txt'), 'a') as f:
                f.write('missing: {}, exception: {}\n'.format(str_id, val))
                continue
        data = base_str_n.loc[str_id].values
        all_n = data[:, np.newaxis] * modifiers[np.newaxis, :]
        outdata[str_id] = _np_describe(all_n)
    outdata = pd.DataFrame(outdata).transpose()
    outdata.to_csv(os.path.join(outdir, 'Ash_Opt_stocastic_n_strs.csv'))


def _np_describe(ndarray, percentiles=(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)):
    """
    pull out the percentiles and put it in a pd.Series
    :param ndarray: a 1 or more dimentional numpy array
    :param percentiles: percentiles to return
    :return: pd.Series
    """
    outdata = pd.Series()
    outdata.loc['mean'] = np.nanmean(ndarray)
    outdata.loc['std'] = np.nanstd(ndarray)
    outdata.loc['min'] = np.nanmin(ndarray)
    for i in percentiles:
        outdata.loc['{}%'.format(int(i * 100))] = np.nanpercentile(ndarray, i * 100)
    outdata.loc['max'] = np.nanmax(ndarray)
    return outdata


def run_all_nload_stuffs():
    """
    wrapper to run all the nloads stocastics
    :return:
    """
    base_outdir = env.gw_met_data(r"mh_modeling\stocastic_n_load_results")

    szdirs = [
        env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely"),
        env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\possible"),
        env.sci(
            r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\probable")
    ]
    with open(env.gw_met_data(r"mh_modeling\stocastic_n_load_results\n_load_names.txt")) as f:
        n_names = []
        for e in f.readlines():
            if e.strip() != '':
                n_names.append(e.strip())
    for sim_end in ['without_trans', 'with_trans']:
        for sz_dir in szdirs:
            for n_name in n_names:
                print('starting N analysis for {} load, {} sims, and {} polygons'.format(n_name, sim_end,
                                                                                         os.path.basename(sz_dir)))
                outdir = os.path.join(base_outdir, sim_end, '{}_{}'.format(n_name, os.path.basename(sz_dir)))
                sims = pd.read_csv(env.gw_met_data(
                    "mh_modeling\stocastic_n_load_results\component_uncertainty_data_{}.csv".format(sim_end)),
                    index_col=0)
                calc_all_ns(sims_org=sims, n_load_name=n_name, outdir=outdir, source_zone_dir=sz_dir)
                output_actual_n_vals(outdir=outdir, mod_dir=outdir)


def spatial_overlays(df1, df2, how='intersection'):  # also in core, but I want a copy in my scripts
    '''Compute overlay intersection of two
        GeoPandasDataFrames df1 and df2
    '''

    df1 = df1.copy()
    df2 = df2.copy()

    if how == 'intersection':
        # Spatial Index to create intersections
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg'] = df1.bbox.apply(lambda x: list(spatial_index.intersection(x)))
        pairs = df1['histreg'].to_dict()
        nei = []
        for i, j in pairs.items():
            for k in j:
                nei.append([i, k])

        pairs = gpd.GeoDataFrame(nei, columns=['idx1', 'idx2'], crs=df1.crs)
        pairs = pairs.merge(df1, left_on='idx1', right_index=True)
        pairs = pairs.merge(df2, left_on='idx2', right_index=True, suffixes=['_1', '_2'])
        pairs['Intersection'] = pairs.apply(lambda x: (x['geometry_1'].intersection(x['geometry_2'])).buffer(0), axis=1)
        pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns, crs=df1.crs)
        cols = pairs.columns.tolist()
        cols.remove('geometry_1')
        cols.remove('geometry_2')
        cols.remove('histreg')
        cols.remove('bbox')
        cols.remove('Intersection')
        dfinter = pairs[cols + ['Intersection']].copy()
        dfinter.rename(columns={'Intersection': 'geometry'}, inplace=True)
        dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=pairs.crs)
        dfinter = dfinter.loc[dfinter.geometry.is_empty == False]
        return (dfinter)
    elif how == 'difference':
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg'] = df1.bbox.apply(lambda x: list(spatial_index.intersection(x)))
        df1['new_g'] = df1.apply(
            lambda x: reduce(lambda x, y: x.difference(y).buffer(0), [x.geometry] + list(df2.iloc[x.histreg].geometry)),
            axis=1)
        df1.geometry = df1.new_g
        df1 = df1.loc[df1.geometry.is_empty == False].copy()
        df1.drop(['bbox', 'histreg', 'new_g'], axis=1, inplace=True)
        return (df1)


if __name__ == '__main__':
    run_all_nload_stuffs()
    print('success, script ran without problems')
