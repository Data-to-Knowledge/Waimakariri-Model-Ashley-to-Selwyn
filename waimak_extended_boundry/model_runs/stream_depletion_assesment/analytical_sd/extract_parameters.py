# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 14/12/2017 11:23 AM
"""

from __future__ import division
import netCDF4 as nc
import numpy as np
from waimak_extended_boundry import smt
from waimak_extended_boundry import \
    get_model
from copy import deepcopy
from waimak_extended_boundry import _get_drn_spd
import pandas as pd
import os


def extract_parameters(outpath, model_id):
    ibnd = smt.get_no_flow()
    m = get_model(model_id)
    outfile = nc.Dataset(outpath, 'w')
    x, y = smt.get_model_x_y(False)
    # create dimensions
    outfile.createDimension('latitude', len(y))
    outfile.createDimension('longitude', len(x))
    outfile.createDimension('layer', smt.layers)

    # create variables
    depth = outfile.createVariable('layer', 'f8', ('layer',), fill_value=np.nan)
    depth.setncatts({'units': 'none',
                     'long_name': 'layer',
                     'missing_value': np.nan})
    depth[:] = range(smt.layers)

    proj = outfile.createVariable('crs', 'i1')  # this works really well...
    proj.setncatts({'grid_mapping_name': "transverse_mercator",
                    'scale_factor_at_central_meridian': 0.9996,
                    'longitude_of_central_meridian': 173.0,
                    'latitude_of_projection_origin': 0.0,
                    'false_easting': 1600000,
                    'false_northing': 10000000,
                    })

    lat = outfile.createVariable('latitude', 'f8', ('latitude',), fill_value=np.nan)
    lat.setncatts({'units': 'NZTM',
                   'long_name': 'latitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_y_coordinate'})
    lat[:] = y

    lon = outfile.createVariable('longitude', 'f8', ('longitude',), fill_value=np.nan)
    lon.setncatts({'units': 'NZTM',
                   'long_name': 'longitude',
                   'missing_value': np.nan,
                   'standard_name': 'projection_x_coordinate'})
    lon[:] = x

    kh = outfile.createVariable('kh', 'f8', ('layer', 'latitude', 'longitude',), fill_value=np.nan)
    kh.setncatts({'units': 'm/day',
                  'long_name': 'horizontal conductivity',
                  'missing_value': np.nan})
    temp = m.upw.hk.array
    temp[ibnd == 0] = np.nan

    kh[:] = temp

    kv = outfile.createVariable('kv', 'f8', ('layer', 'latitude', 'longitude',), fill_value=np.nan)
    kv.setncatts({'units': 'm/day',
                  'long_name': 'vertical conductivity',
                  'missing_value': np.nan})
    temp = m.upw.vka.array
    temp[ibnd == 0] = np.nan
    kv[:] = temp

    cond = outfile.createVariable('str_cond', 'f8', ('latitude', 'longitude',), fill_value=np.nan)
    cond.setncatts({'units': 'm2/day',
                    'long_name': 'stream/drain conductance',
                    'missing_value': np.nan})

    drn = _get_drn_spd(smt.reach_v, smt.wel_version, n_car_dns=False)
    temp = smt.df_to_array(drn, 'i')
    idx = np.isfinite(temp)
    drn_cond = m.drn.stress_period_data.array['cond'][0, 0]
    drn_width = deepcopy(drn_cond)
    drn_width[np.isfinite(drn_width)] = 200
    drn_len = deepcopy(drn_width)
    drn_cond[~idx] = np.nan
    drn_width[~idx] = np.nan
    drn_len[~idx] = np.nan

    sfr_len = m.sfr.stress_period_data.array['rchlen'][0, 0]
    sfr_seg = m.sfr.stress_period_data.array['iseg'][0, 0]
    temp = pd.DataFrame(m.sfr.segment_data[0]).set_index('nseg').loc[:, 'width1'].to_dict()
    temp[np.nan] = np.nan
    sfr_width = smt.get_empty_model_grid() * np.nan
    sfr_width[np.isfinite(sfr_seg)] = vec_translate(sfr_seg[np.isfinite(sfr_seg)], temp)
    sfr_k = m.sfr.stress_period_data.array['strhc1'][0, 0]
    sfr_cond = sfr_k * sfr_width * sfr_len
    drn_cond[np.isfinite(sfr_cond)] = sfr_cond[np.isfinite(sfr_cond)]
    cond[:] = drn_cond


def extract_parameters_to_points(outdir, model_id):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ibnd = smt.get_no_flow()
    m = get_model(model_id)
    x, y = smt.get_model_x_y()
    layers = range(smt.layers)
    elv_db = smt.calc_elv_db()
    thickness = elv_db[:-1] - elv_db[1:]

    kh = m.upw.hk.array
    kh[ibnd == 0] = np.nan

    kv = m.upw.vka.array
    kv[ibnd == 0] = np.nan

    for l in layers:
        idx = np.isfinite(kv[l])
        temp = pd.DataFrame({'x': x[idx].flatten(), 'y': y[idx].flatten(),
                             'kh': kh[l][idx].flatten(), 'kv': kv[l][idx].flatten(),
                             'thickness': thickness[l][idx].flatten()})
        temp.to_csv(os.path.join(outdir, 'kh_kv_thick_layer_{:02d}.csv'.format(l + 1)))

    drn = _get_drn_spd(smt.reach_v, smt.wel_version, n_car_dns=False)
    temp = smt.df_to_array(drn, 'i')
    idx = np.isfinite(temp)
    drn_cond = m.drn.stress_period_data.array['cond'][0, 0]
    drn_width = deepcopy(drn_cond)
    drn_width[np.isfinite(drn_width)] = 200
    drn_len = deepcopy(drn_width)
    drn_cond[~idx] = np.nan
    drn_width[~idx] = np.nan
    drn_len[~idx] = np.nan

    sfr_len = m.sfr.stress_period_data.array['rchlen'][0, 0]
    sfr_seg = m.sfr.stress_period_data.array['iseg'][0, 0]
    temp = pd.DataFrame(m.sfr.segment_data[0]).set_index('nseg').loc[:, 'width1'].to_dict()
    temp[np.nan] = np.nan
    sfr_width = smt.get_empty_model_grid() * np.nan
    sfr_width[np.isfinite(sfr_seg)] = vec_translate(sfr_seg[np.isfinite(sfr_seg)], temp)
    sfr_k = m.sfr.stress_period_data.array['strhc1'][0, 0]
    sfr_cond = sfr_k * sfr_width * sfr_len
    drn_cond[np.isfinite(sfr_cond)] = sfr_cond[np.isfinite(sfr_cond)]

    idx = np.isfinite(drn_cond)
    temp = pd.DataFrame({'x': x[idx].flatten(), 'y': y[idx].flatten(),
                         'cond':drn_cond[idx].flatten()})
    temp.to_csv(os.path.join(outdir, 'conductance.csv'))


def vec_translate(a, d):
    return np.vectorize(d.__getitem__)(a)


if __name__ == '__main__':
    extract_parameters_to_points(
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\stream_depletion\point_parameters",
        'NsmcBase')
