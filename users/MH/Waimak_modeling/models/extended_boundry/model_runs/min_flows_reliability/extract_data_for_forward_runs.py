# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 28/09/2017 2:00 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction import \
    get_hds_at_wells, get_flow_at_points, get_samp_points_df
from glob import glob
import os
import datetime
from copy import deepcopy
import flopy
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.convergance_check import \
    zipped_modflow_converged
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
from visualise_data_from_fruns import gw_site_groups
import itertools
import netCDF4 as nc
import sys


def extract_forward_run(name_file_path):
    """
    extract the data from a forward run
    :param name_file_path: path to teh mf name file with or without extension
    :return: pd.DataFrame
    """
    wells = list(itertools.chain(*gw_site_groups.values())) # set in visualise data
    streams = get_samp_points_df()
    streams = streams.loc[streams.m_type == 'min_flow'].index

    hd_data = get_hds_at_wells(wells, rel_kstpkpers=0, name_file_path=name_file_path,
                               add_loc=True)
    hd_data = hd_data.rename(columns={'hd_m3d_kstp0_kper0': 'kstpkper0'})
    str_data = get_flow_at_points(streams, rel_kstpkpers=0, base_path=name_file_path)
    str_data = str_data.rename(columns={'flow_m3d_kstp0_kper0': 'kstpkper0'})
    outdata = pd.concat((hd_data, str_data))
    return outdata


def extract_forward_metadata(forward_run_dir, outpath):
    """
    extract and save the metadata senario, convergence, path, name, rcp ect...
    :param forward_run_dir:
    :param outpath: path to save the data at
    :return: outpath
    """
    paths = glob(os.path.join(forward_run_dir, '*/*.nam'))
    model_id = os.path.basename(paths[0]).split('_')[0]
    outpath = os.path.join(os.path.dirname(outpath), '{}_{}'.format(model_id, os.path.basename(outpath)))
    model_names = [os.path.basename(path).replace('.nam', '').replace(model_id + '_', '') for path in paths]
    converged = [zipped_modflow_converged(path, return_nans=True) for path in paths]
    outdata = pd.DataFrame(index=model_names, data={'converged': converged, 'path': paths})

    # extract the sen at ect from the path?
    outdata['is_cc'] = idx = outdata.index.str.contains('RCP')
    outdata.loc[outdata.index.str.contains('low_3_m'), 'amalg_type'] = 'low_3_m'
    outdata.loc[outdata.index.str.contains('tym'), 'amalg_type'] = 'tym'
    outdata.loc[idx, 'rcm'] = [e.replace('_' + outdata.loc[e, 'amalg_type'], '').split('_')[-3] for e in
                               outdata.index[idx]]
    outdata.loc[idx, 'rcp'] = [e.replace('_' + outdata.loc[e, 'amalg_type'], '').split('_')[-2] for e in
                               outdata.index[idx]]
    outdata.loc[idx, 'period'] = [e.replace('_' + outdata.loc[e, 'amalg_type'], '').split('_')[-1] for e in
                                  outdata.index[idx]]
    outdata.loc[idx, 'sen'] = [e.replace(
        '_{rcm}_{rcp}_{period}_{amalg_type}'.format(**outdata.loc[e, ['rcm', 'rcp', 'period', 'amalg_type']].to_dict()),
        '') for e in outdata.index[idx]]
    outdata.loc[~idx, 'sen'] = outdata.loc[~idx].index

    outdata.to_csv(outpath)
    return outpath

def extract_and_save_all_cc_mult_missing_w(forward_run_dir, outpath):
    """
    extract and save all of the ccmults and missing waters from teh forward runs
    :param forward_run_dir: directory with teh forward runs
    :param outpath: path to save csv at
    :return: outpath
    """
    paths = glob(os.path.join(forward_run_dir, '*/cc_mult_info.txt'))
    paths = sorted(paths)
    model_id = os.path.basename(os.path.dirname(paths[0])).split('_')[0]
    outpath = os.path.join(os.path.dirname(outpath), '{}_{}'.format(model_id, os.path.basename(outpath)))
    temp_names = [os.path.basename(os.path.dirname(e)).replace(model_id + '_', '') for e in paths]
    outdata = pd.DataFrame(index=temp_names,columns=['cc_mult','missing_water'])
    for name, path in zip(temp_names, paths):
        with open(path) as f:
            lines = f.readlines()
            if len(lines) != 2:
                raise ValueError('more than 2 lines in {}'.format(path))
            outdata.loc[name,'missing_water'] = float(lines[1].split(':')[-1].strip())
            outdata.loc[name, 'cc_mult'] = float(lines[0].split(':')[-1].strip())

    with open(outpath, 'w') as f:
        f.write('cc_mult and missing water from model {}. cc_mult is unitless and is applied to the irrigation '
                'pumping wells in the waimakariri zone.  pumping is truncated at the max CAV and allocation. '
                'missing water is the increase in irrgation demand for the LSR layer that cannot be met through '
                'groundwater and is in units of m3/day; made {}\n'.format(model_id,
                                                                          datetime.datetime.now().isoformat()))

    outdata.to_csv(outpath, mode='a')
    return outpath



def extract_and_save_all_forward_runs(forward_run_dir, outpath):
    """
    extract and save all of the forward run data as absolute flow and hds
    :param forward_run_dir: directory with the forward run
    :param outpath: the path to save the csv
    :return: outpath
    """
    paths = glob(os.path.join(forward_run_dir, '*/*.nam'))
    paths = sorted(paths)
    model_id = os.path.basename(paths[0]).split('_')[0]
    outpath = os.path.join(os.path.dirname(outpath), '{}_{}'.format(model_id, os.path.basename(outpath)))
    for i, path in enumerate(paths):
        temp_name = os.path.basename(path).replace('.nam', '').replace(model_id + '_', '')
        temp = extract_forward_run(path)
        temp = temp.rename(columns={'kstpkper0': temp_name})
        if i == 0:
            outdata = temp.loc[:, ['depth', 'i', 'j', 'k', 'mid_screen_elv', 'nztmx', 'nztmy', temp_name]]
        else:
            temp = temp.drop(['depth', 'i', 'j', 'k', 'mid_screen_elv', 'nztmx', 'nztmy'], 1)
            outdata = pd.merge(outdata, temp, right_index=True, left_index=True)
    with open(outpath, 'w') as f:
        f.write('flow and flux values from model {}. all flow values in m3/day; all hd; z in m; x; y in nztm, '
                'i;j;k are unit less; made {}\n'.format(model_id, datetime.datetime.now().isoformat()))
    outdata.index.name = 'site'
    outdata.to_csv(outpath, mode='a')
    return outpath


def make_rel_data(data_path, meta_data_path, out_path):
    """
    convert the absolute data to reletive data
    :param data_path: absolute data path
    :param meta_data_path: path to the meta data
    :param out_path: path to save the relative data
    :return:
    """
    org_data = pd.read_csv(data_path, skiprows=1, index_col=0)
    meta_data = pd.read_csv(meta_data_path, index_col=0)
    model_id = os.path.basename(data_path).split('_')[0]
    out_path = os.path.join(os.path.dirname(out_path), '{}_{}'.format(model_id, os.path.basename(out_path)))
    remove_keys = ['depth', 'i', 'j', 'k', 'mid_screen_elv', 'nztmx', 'nztmy']
    divide_keys = list(set(org_data.keys()) - set(remove_keys))
    outdata = deepcopy(org_data)
    actual_keys = []
    well_idx = outdata.index.str.contains('/')
    for key in divide_keys:
        reference_key = get_baseline_name(meta_data=meta_data, name=key)

        if reference_key == key:  # do not make reference keys relative
            actual_keys.append(key)
            continue

        # set up indexes for dry/no_flow in sim/ref
        dry_ref = ((np.isclose(org_data.loc[:, reference_key], -888) & well_idx) |
                   (np.isclose(org_data.loc[:, reference_key], 0) & ~well_idx))
        dry_sim = (np.isclose(outdata.loc[:, key], -888) & (well_idx))

        # relative streams
        outdata.loc[~well_idx, key] *= 1 / org_data.loc[~well_idx, reference_key]
        # actual drawdown
        outdata.loc[well_idx, key] = outdata.loc[well_idx, key] - org_data.loc[well_idx, reference_key]

        # keys for dry well -888, dry(well) or zero flow (sfr/drn) in reference (-777) zero flow is left as 0%
        outdata.loc[dry_sim, key] = -888
        outdata.loc[dry_ref, key] = -777

    # write data with a header
    with open(out_path, 'w') as f:
        f.write(
            'relative flow and flux values from model {mid}. {act} are actual data'
             'all other values relative: '
             'any climate change senario is relative to the mean of RCPpast for the senario and rcm'
             'all others are relative to current'
             'draw down (wells) is senario - baseline'
             'for actuals all flow values in m3/day; all hd; z in m; x; y in nztm'
             'i;j;k are unit less; dry wells are filled with a value of -888 and '
            'dry (well) or noflow (stream) in reference simulation is -777 '
            'made: {dt}\n'.format(mid=model_id, act=str(actual_keys).replace(',',';'),
                                                         dt=datetime.datetime.now().isoformat()))
    # write data
    outdata.to_csv(out_path, mode='a')


def get_baseline_name(meta_data, name, raise_non_converged=True):
    """
    get teh appropriate baseline to make relative data
    cc_senarios: rcppast for the senario and rcm
    non_cc_senarios: the current (2014/15 pumping senario)
    :param meta_data: the metadata pd.DataFrame
    :param name: the name of the run
    :param raise_non_converged: if True raise a value error if the baseline senario did not converge
    :return:
    """
    if not meta_data.loc[name, 'is_cc']:
        outname = 'current'

    else:
        rcp = 'RCPpast'
        sen = meta_data.loc[name, 'sen']
        amalg_type = 'tym'
        rcm = meta_data.loc[name, 'rcm']
        period = 1980
        outname= '{}_{}_{}_{}_{}'.format(sen, rcm, rcp, period, amalg_type)

    if raise_non_converged:
        if pd.isnull(meta_data.loc[outname, 'converged']):
            raise ValueError('null convergence for reference scenario {}'.format(outname))
        if not meta_data.loc[outname, 'converged']:
            raise ValueError('reference scenario {} did not converge'.format(outname))

    return outname


def plt_drawdown(meta_data_path, outdir,raise_non_converged=True):
    """
    plot the drawdown for each senario as compared to the base senario
    :param meta_data_path: path to the metadata
    :param outdir: directory to save the plots
    :param raise_non_converged: passed to get baseline senario
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    meta_data = pd.read_csv(meta_data_path, index_col=0)
    for name in meta_data.index:
        converged = meta_data.loc[name, 'converged']
        if pd.isnull(converged):
            continue
        ref_name = get_baseline_name(meta_data, name, raise_non_converged)
        mod_per_hds_path =meta_data.loc[ref_name, 'path'].replace('.nam', '.hds')
        mod_per_hds = flopy.utils.HeadFile(mod_per_hds_path).get_data((0, 0))
        plt_out_dir = os.path.join(outdir, name)
        if not os.path.exists(plt_out_dir):
            os.makedirs(plt_out_dir)

        hd_file_path = meta_data.loc[name, 'path'].replace('.nam', '.hds')
        hds = flopy.utils.HeadFile(hd_file_path).get_data((0, 0))
        # set no flow and dry cells to nan, in plotting, the drycells with appear green, the others will appear black
        hds[hds > 1e20] = np.nan
        hds[hds < -666] = np.nan
        hds = hds - mod_per_hds
        for layer in range(smt.layers):
            fig, ax = smt.plt_matrix(hds[layer], vmin=-5, vmax=5,
                                     title='draw_down for {}, converged: {}, layer {:02d}'.format(name, converged,
                                                                                                  layer),
                                     no_flow_layer=layer, cmap='RdBu')
            ax.set_facecolor('g')
            fig.savefig(os.path.join(plt_out_dir, 'drawdown_layer_{:02d}.png'.format(layer)))
            plt.close(fig)

def netcdf_drawdown (meta_data_path, outpath, readme, raise_non_converged=True):

    meta_data = pd.read_csv(meta_data_path, index_col=0)
    data = {}
    for name in meta_data.index:
        converged = meta_data.loc[name, 'converged']
        if pd.isnull(converged):
            continue
        ref_name = get_baseline_name(meta_data, name, raise_non_converged)
        mod_per_hds_path =meta_data.loc[ref_name, 'path'].replace('.nam', '.hds')
        mod_per_hds = flopy.utils.HeadFile(mod_per_hds_path).get_data((0, 0))
        mod_per_hds[mod_per_hds > 1e20] = np.nan
        mod_per_hds[mod_per_hds < -666] = np.nan

        hd_file_path = meta_data.loc[name, 'path'].replace('.nam', '.hds')
        hds = flopy.utils.HeadFile(hd_file_path).get_data((0, 0))
        # set no flow and dry cells to nan, in plotting, the drycells with appear green, the others will appear black
        hds[hds > 1e20] = np.nan
        hds[hds < -666] = np.nan
        hds = hds - mod_per_hds # drawdown appears negative
        data[name] = hds

    # for each loop create a netcdf file
    outfile = nc.Dataset(outpath, 'w')
    x, y = smt.get_model_x_y(False)
    # create dimensions
    outfile.createDimension('latitude', len(y))
    outfile.createDimension('longitude', len(x))
    outfile.createDimension('layer', smt.layers - 1)

    # create variables
    depth = outfile.createVariable('layer', 'f8', ('layer',), fill_value=np.nan)
    depth.setncatts({'units': 'none',
                     'long_name': 'layer',
                     'missing_value': np.nan,
                     'comments': '0 indexed'})
    depth[:] = range(smt.layers)

    proj = outfile.createVariable('crs', 'i1')
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

    for name, values in data.items():
        temp = outfile.createVariable(name, 'f8', ('layer','latitude','longitude',),
                                      fill_value=np.nan)
        temp.setncatts({'units': 'm',
                       'long_name': name,
                       'missing_value': np.nan})
        temp[:] = values

    outfile.description = (
        'drawdown (negitive is drawdown vs reference) for forward runs: {}'.format(readme))
    outfile.history = 'created {}'.format(datetime.datetime.now().isoformat())
    outfile.source = 'script: {}'.format(sys.argv[0])
    outfile.close()
    #todo check


def gen_all_outdata_forward_runs(forward_run_dir, outdir, plt_dd=False):
    """
    a wrapper to extract all data from the forward runs
    :param forward_run_dir: the directory that holds the forward runs
    :param outdir: the directory to save all teh data in
    :param plt_dd: boolean if True plot the drawdown for each senario
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    absolute_outpath = 'absolute_data.csv'
    meta_data_path = 'meta_data.csv'
    relative_outpath = 'relative_data.csv'
    cc_mult_path = 'cc_mult_miss_water.csv'
    netcdf_path = 'drawdown.nc'

    # copy readme over
    with open(os.path.join(forward_run_dir,'READ_ME.txt')) as f:
        readme = f.read()
    with open(os.path.join(outdir,'READ_ME.txt'),'w') as f:
        f.write(readme)

    print('extracting cc_mult and missing water')
    extract_and_save_all_cc_mult_missing_w(forward_run_dir,os.path.join(outdir, cc_mult_path))

    print('extracting absolute data')
    absolute_outpath = extract_and_save_all_forward_runs(forward_run_dir, os.path.join(outdir, absolute_outpath))

    print('extracting metadata')
    meta_data_path = extract_forward_metadata(forward_run_dir, os.path.join(outdir, meta_data_path))

    print('creating relative data')
    make_rel_data(absolute_outpath, meta_data_path, os.path.join(outdir, relative_outpath))

    print('creating netcdf of drawdown')
    netcdf_drawdown(meta_data_path, os.path.join(outdir,netcdf_path),readme)

    if plt_dd:
        print('plotting drawdown')
        plt_drawdown(meta_data_path, os.path.join(outdir, 'plots'))


if __name__ == '__main__':
    #tests (run in the run script for forward runs)
    gen_all_outdata_forward_runs(
        r"D:\mh_model_runs\forward_runs_2017_10_10",
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\forward_sw_gw\results\cc_only_to_waimak",
        True)
    extract_and_save_all_cc_mult_missing_w(r"C:\Users\MattH\Desktop\forward_run_test",r"C:\Users\MattH\Downloads\test_ccmult_extract.csv")
    print('done')
