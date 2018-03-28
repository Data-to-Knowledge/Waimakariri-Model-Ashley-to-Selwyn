# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/03/2018 3:15 PM
"""

from __future__ import division
from core import env
from hydraulics_cleaned import hunt2003, theis_jenkins
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import itertools
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col
from glob import glob
import os
import geopandas as gpd
import flopy
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.stream_depletion_assesment.stream_depletion_numerical_model_runs.starting_hds_ss_sy import \
    get_sd_well_list
import matplotlib.pyplot as plt


# is it worth doing a bulk comparison for each stream (e.g. all reaches?)
# I should consider dry reaches as well
# point to point nodes for all shapefiles might be the simplest
# consider an assumed pumped aquifer thickness of 25m starting from the well midpoint.  area above is the 'aquitard'

def get_transmissivity_kaqt(row, col, layer, pumping_thickness, khs, kvs, radius):
    num_cells = radius // smt.grid_space
    rows = np.arange(row - num_cells, row + num_cells + 1).astype(int)
    rows = rows[(rows >= 0) & (rows < smt.rows)]  # get rid of all rows, and cols which are out of order
    cols = np.arange(col - num_cells, col + num_cells + 1).astype(int)
    cols = cols[(cols >= 0) & (cols < smt.cols)]
    idxs = np.array(list(itertools.product([int(layer)], rows, cols)))
    idxs = [idxs[:, 0], idxs[:, 1], idxs[:, 2]]
    trans = np.nanmean(khs[idxs]) * pumping_thickness

    layers = [0]
    layers.extend(range(1, int(layer)))
    idxs = np.array(list(itertools.product(layers, rows, cols)))
    idxs = [idxs[:, 0], idxs[:, 1], idxs[:, 2]]
    kaqt = np.nanmean(kvs[idxs])
    return trans, kaqt


def get_baqt(layer, row, col, mid_screen_elv, elv_db, heads):
    layer = int(layer)
    row = int(row)
    col = int(col)
    cell_bot = elv_db[layer + 1, row, col]
    sat_thickness = np.nan
    if heads[layer, row, col] < cell_bot:  # cell is dry what then
        sat_thickness = abs((elv_db[0, row, col] - mid_screen_elv) / 2)
        # abs in the unlikely chance that the well is slightly above ground level #todo talk over
    else:
        idx = heads[:, row, col] > elv_db[1:, row, col]
        if not idx.any():
            sat_thickness = np.nan
        else:
            top_sat_layer = (np.arange(11)[idx]).min()
            sat_thickness = min(elv_db[top_sat_layer, row, col], heads[top_sat_layer, row, col]) - mid_screen_elv

    if np.isnan(sat_thickness):
        pass
    elif sat_thickness <= 0:
        sat_thickness = abs((elv_db[0, row, col] - mid_screen_elv) / 2)
        # abs in the unlikely chance that the well is slightly above ground level #todo talk over

    return sat_thickness


def calc_analytical_sd(well_nums, name_file_path, outpath, radius=5000):
    """
    calculate the hunt(2003) stream depletion for every well and swaz.
    T is defined by an assumed 25m thick aquifer and
    :param well_nums:
    :param name_file_path:
    :param outpath:
    :param radius:
    :return:
    """
    no_flow = smt.get_no_flow()
    m = flopy.modflow.Modflow.load(name_file_path, model_ws=os.path.dirname(name_file_path), forgive=False)
    # get well numbers/location
    all_wells = get_all_well_row_col()
    data = all_wells.loc[well_nums, [u'nztmx', u'nztmy', u'depth', u'mid_screen_depth',
                                     u'mid_screen_elv', u'row', u'col', u'layer']]
    data = data.dropna()
    # i do not need to actually check dry reaches as they are inherent in the swazes.

    # get times
    times = [7, 150]
    s = 0.0001  # todo assumed pumped s
    sy = 0.1  # todo assumed aquitard sy

    # get transmissivity, kaqt
    kvs = m.upw.vka.array.copy()
    kvs[no_flow == 0] = np.nan
    khs = m.upw.hk.array.copy()
    khs[no_flow == 0] = np.nan
    pumping_thickness = 25  # assume 25 m pumped aquifer
    trans, kaqt = np.vectorize(get_transmissivity_kaqt,
                               excluded={'khs', 'kvs', 'radius', 'pumping_thickness'},
                               otypes=[np.float, np.float])(row=data.row, col=data.col,
                                                            layer=data.layer,
                                                            pumping_thickness=pumping_thickness,
                                                            kvs=kvs,
                                                            khs=khs, radius=radius)
    data.loc[:, 'trans'] = trans
    data.loc[:, 'kaqt'] = kaqt

    # get baqt assume from midscreen elevation to top of saturated zone
    elv_db = smt.calc_elv_db()
    heads = flopy.utils.HeadFile(name_file_path.replace('.nam', '.hds')).get_alldata()[0]
    heads[heads > 1e15] = np.nan
    data.loc[:, 'baqt'] = np.vectorize(get_baqt, excluded={'elv_db', 'heads'},
                                       otypes=[np.float])(layer=data.layer,
                                                          row=data.row,
                                                          col=data.col,
                                                          mid_screen_elv=data.mid_screen_elv,
                                                          elv_db=elv_db,
                                                          heads=heads)

    # get separation distance
    # stream conductance for each name
    swaz_paths = glob(os.path.join(smt.sdp, r"m_ex_bd_inputs\raw_sw_samp_points\*\swaz\*.shp"))
    swaz_names = [os.path.basename(sp).replace('.shp', '') for sp in swaz_paths]
    stream_conducts = {}
    all_str_conducts = m.drn.stress_period_data.array['cond'][0, 0]
    seg_to_width = dict(zip(m.sfr.segment_data[0]['nseg'], m.sfr.segment_data[0]['width1']))
    temp = m.sfr.stress_period_data.array
    seg_to_width[-99] = np.nan
    seg = temp['iseg'][0, 0]
    seg[np.isnan(seg)] = -99
    temp2 = temp['strhc1'][0, 0] * temp['rchlen'][0, 0] * np.vectorize(seg_to_width.__getitem__)(seg)
    all_str_conducts[np.isfinite(temp2)] = temp2[np.isfinite(temp2)]
    for sp, name in zip(swaz_paths, swaz_names):
        geo = gpd.read_file(sp)
        geo_xs = geo.geometry.x[np.newaxis, :]
        geo_ys = geo.geometry.y[np.newaxis, :]
        xs = data.loc[:, 'nztmx'].values[:, np.newaxis] - geo_xs
        ys = data.loc[:, 'nztmy'].values[:, np.newaxis] - geo_ys  # todo check size
        data.loc[:, 'sep_dist_{}'.format(name)] = ((xs ** 2 + ys ** 2) ** 0.5).min(axis=1)
        idx = np.isfinite(smt.shape_file_to_model_array(sp, 'k', True))
        stream_conducts[name] = np.nanmean(all_str_conducts[idx])

    for time in times:
        for name in swaz_names:
            sd = hunt2003(discharge=100,
                          # discharge is not relevant as I am returning the fraction rather than the rate
                          time=time,
                          trans=data.loc[:, 'trans'],
                          s=s,
                          kaqt=data.loc[:, 'kaqt'],
                          baqt=data.loc[:, 'baqt'],
                          sy=sy,
                          lam=stream_conducts[name],
                          l=data.loc[:, 'sep_dist_{}'.format(name)],
                          return_rate=False)
            data.loc[:, 'hunt_{}_sd_{}'.format(name, time)] = sd

            tsd = theis_jenkins(discharge=100,
                                time = time,
                                trans=data.loc[:,'trans'],
                                s=sy,
                                l=data.loc[:, 'sep_dist_{}'.format(name)],
                                return_rate=False)
            data.loc[:, 'theis_{}_sd_{}'.format(name, time)] = tsd

    if not os.path.exists(os.path.dirname(outpath)):
        os.makedirs(os.path.dirname(outpath))
    data.index.name = 'well'

    data.loc[:,'nearest_swaz'] = np.vectorize(_get_closest_stream,
                                              excluded=['streams'])(streams=np.array(swaz_names),
                                                                    distances=[_oblist(e) for e in data.loc[:,['sep_dist_{}'.format(name) for name in swaz_names]].values])
    for time in times:
        for method in ['hunt', 'theis']:
            for well in data.index:
                stream = data.loc[well, 'nearest_swaz']
                data.loc[well,'{}_nearest_sd_{}'.format(method,time)] = data.loc[well, '{}_{}_sd_{}'.format(method, stream, time)]

    data.to_csv(outpath)
    return data

def _get_closest_stream(streams, distances):
    idx = np.argmin(distances.data),
    return streams[idx]

class _oblist(object):
    def __init__(self,l):
        self.data = l

def join_sds(numerical_sd_7, numerical_sd_150, analytical_data_hunt):
    # make the data into long mode
    an_keys = analytical_data_hunt.keys()
    an_keys = an_keys[an_keys.str.contains('sd_7')] #todo update for theis addition
    sd7_an = pd.melt(analytical_data_hunt.loc[:, an_keys].reset_index(), id_vars='well', var_name='stream',
                     value_name='analytical_hunt')
    sd7_an.loc[:, 'stream'] = sd7_an.loc[:, 'stream'].str.replace('_sd_7', '')
    sd7_an.loc[:, 'analytical_hunt'] *=100

    num_keys = ['custmaindrain_swaz',
                'cust_swaz',
                'eyre_swaz',
                'n7drain_swaz',
                'cam_swaz',
                'courtenay_swaz',
                'greigs_swaz',
                'kaiapoi_swaz',
                'kairaki_swaz',
                'northbrook_swaz',
                'ohoka_swaz',
                'saltwater_swaz',
                'southbrook_swaz',
                'taranaki_swaz',
                'waikuku_swaz',
                'ashley_swaz',
                'waimakupper_swaz',
                'waimaklower_swaz',
                'waimak_swaz'
                ]
    sd7_num = pd.melt(numerical_sd_7.loc[:, num_keys].reset_index(), id_vars='well', var_name='stream',
                      value_name='numerical')

    sd7 = pd.merge(sd7_an, sd7_num, left_on=['well', 'stream'], right_on=['well', 'stream'])

    # sd150
    an_keys = analytical_data_hunt.keys()
    an_keys = an_keys[an_keys.str.contains('sd_150')]
    sd150_an = pd.melt(analytical_data_hunt.loc[:, an_keys].reset_index(), id_vars='well', var_name='stream',
                       value_name='analytical_hunt')
    sd150_an.loc[:, 'stream'] = sd150_an.loc[:, 'stream'].str.replace('_sd_150','')

    sd150_num = pd.melt(numerical_sd_150.loc[:, num_keys].reset_index(), id_vars='well', var_name='stream',
                        value_name='numerical')
    sd150_an.loc[:, 'analytical_hunt'] *=100

    sd150 = pd.merge(sd150_an, sd150_num, left_on=['well', 'stream'], right_on=['well', 'stream'])
    return sd7, sd150

# todo add a comparison for my hks and matt's khs


if __name__ == '__main__':
    well_nums = get_sd_well_list('NsmcBase')
    analytical_data = calc_analytical_sd(well_nums,
                       r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\NsmcBase\AW20171024_2_i2_optver\i2\mf_aw_ex.nam",
                       r"C:\Users\MattH\Downloads\analytical_sd.csv")
    numerical_sd_7 = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\med_s\NsmcBase_2017_12_26_data\NsmcBase_extract_sd7.csv", index_col=0, skiprows=1)
    numerical_sd_150 = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\numerical_results\med_s\NsmcBase_2017_12_26_data\NsmcBase_extract_sd150.csv", index_col=0, skiprows=1)
    sd7, sd150 = join_sds(numerical_sd_7, numerical_sd_150, analytical_data)
    plt.scatter(sd7.analytical, sd7.numerical, c='red', label='sd7')
    plt.scatter(sd150.analytical, sd150.numerical, c='blue', label='sd150')
    plt.legend()
    plt.xlabel('analytical')
    plt.ylabel('numerical')
    plt.show()
    print('done')

