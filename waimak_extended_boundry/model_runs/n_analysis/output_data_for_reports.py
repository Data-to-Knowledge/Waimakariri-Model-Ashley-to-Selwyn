# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/05/2018 8:57 AM
"""

from __future__ import division
import env
import os
import pandas as pd
import numpy as np
import itertools
from get_current_pathway_n import get_pc5pa_mults
from waimak_extended_boundry import smt
from percentage_reduction_maps import wdc_wells, private_wells, streams, gen_stream_targets, gen_well_targets
from glob import glob


def output_current_pathways_table(data_path, outdir, prefix=''):  # todo
    data = pd.read_csv(data_path, index_col=0, header=[0, 1])

    # see n waimak  but basically CMP, GMP, PC5PA, current pathways (half PC5PA)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    current_measured = gen_well_targets('current_measured')
    current_measured.update(gen_stream_targets('current_measured'))
    temp_names = ['current', 'cmp', 'gmp', 'pc5pa', 'cpath']
    pc5_mults = get_pc5pa_mults()

    for recptor_name in ['wdc_wells', 'private_wells', 'streams']:
        receptor_set = sorted(eval(recptor_name))
        outpath = os.path.join(outdir, '{}current_pathways_{}.csv'.format(prefix, recptor_name))
        with open(outpath, 'w') as f:
            f.write('Site,Current_measured (mg/L),CMP_N (mg/L),GMP_N (mg/L),PC5PA_N (mg/L),Current_pathways_N (mg/L)\n')
        outdata = pd.DataFrame(index=receptor_set, columns=temp_names)
        outdata.loc[:, 'current'] = [current_measured[e] for e in outdata.index]

        # cmp
        temp = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for
                f, m, n in data['cmp'].loc[receptor_set, ['5%', '50%', '95%']].itertuples(False, None)]
        outdata.loc[:, 'cmp'] = temp

        # gmp
        gmp = data['gmp'].loc[receptor_set, ['5%', '50%', '95%']]
        temp = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for f, m, n in gmp.itertuples(False, None)]
        outdata.loc[:, 'gmp'] = temp

        # pc5pa
        pc5 = (gmp.transpose() * (1 + pc5_mults.loc[gmp.index])).transpose()
        temp = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for f, m, n in pc5.itertuples(False, None)]
        outdata.loc[:, 'pc5pa'] = temp

        # cpath
        cpath = (gmp.transpose() * (1 + pc5_mults.loc[gmp.index] / 2)).transpose()
        temp = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for f, m, n in cpath.itertuples(False, None)]
        outdata.loc[:, 'cpath'] = temp

        outdata.to_csv(outpath, header=False, mode='a')


def output_future_scenarios_table_withpa(scenarios, data_path, outdir, prefix=''):
    data = pd.read_csv(data_path, index_col=[0], header=[0, 1])
    if isinstance(scenarios, str):
        if scenarios == 'all':
            scenarios = list(data.keys().levels[0])
    assert np.in1d(scenarios, data.keys().levels[0]).all(), 'missing scenarios in data file'
    pc5_mults = get_pc5pa_mults()

    for recptor_name in ['wdc_wells', 'private_wells', 'streams']:
        receptor_set = sorted(eval(recptor_name))
        outpath = os.path.join(outdir, '{}projection_with_pa_{}.csv'.format(prefix, recptor_name))
        outdata = pd.DataFrame(index=receptor_set, columns=scenarios)

        for scen in scenarios:
            temp = data[scen].loc[receptor_set, ['5%', '50%', '95%']]
            temp = (temp.transpose() * (1 + pc5_mults.loc[temp.index] / 2)).transpose()
            temp2 = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for f, m, n in
                     temp.itertuples(False, None)]
            outdata.loc[:, scen] = temp2
        outdata.to_csv(outpath)


def output_future_scenarios_table_without(scenarios, data_path, outdir, prefix=''):
    data = pd.read_csv(data_path, index_col=[0], header=[0, 1])
    if isinstance(scenarios, str):
        if scenarios == 'all':
            scenarios = list(data.keys().levels[0])
    assert np.in1d(scenarios, data.keys().levels[0]).all(), 'missing scenarios in data file'

    for recptor_name in ['wdc_wells', 'private_wells', 'streams']:
        receptor_set = sorted(eval(recptor_name))
        outpath = os.path.join(outdir, '{}projection_without_pa_{}.csv'.format(prefix, recptor_name))
        outdata = pd.DataFrame(index=receptor_set, columns=scenarios)

        for scen in scenarios:
            temp = data[scen].loc[receptor_set, ['5%', '50%', '95%']]
            temp2 = ['{} ({}-{})'.format(round(m, 1), round(f, 1), round(n, 1)) for f, m, n in
                     temp.itertuples(False, None)]
            outdata.loc[:, scen] = temp2
        outdata.to_csv(outpath)


def output_interzone_table(scenarios, data_path, outpath):
    """
    CREATE INTERZONE TABLE AND ADD CHCH 8 TO THE DATA
    :param scenarios:
    :param data_path:
    :param outpath:
    :return:
    """
    data = pd.read_csv(data_path, index_col=[0, 1], header=[0, 1])
    if isinstance(scenarios, str):
        if scenarios == 'all':
            scenarios = list(data.keys().levels[0])
    assert np.in1d(scenarios, data.keys().levels[0]).all(), 'missing scenarios in data file'

    chch_8 = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and re"
                                 r"sults\ex_bd_va\n_results\interzone_n_results\n_data.csv"),
                         index_col=[0, 1], header=[0, 1])['chch_8kgha']
    percentiles = ['5%', '50%', '95%', '99%']
    out_keys = 'statistics, 5th percentile, Median, 95th percentile, 99th percentile\n'
    with open(outpath, 'w') as f:
        f.write(out_keys)
        f.write('Shallow aquifer nitrate-N (mg/L)\n')

    # shallow (riccaton)
    outdata = pd.DataFrame(index=scenarios, columns=percentiles)
    for scen, per in itertools.product(scenarios, percentiles):
        outdata.loc[scen, per] = data.loc[('Riccaton', 'full_city'), (scen, per)] + chch_8.loc[
            ('Riccaton', 'full_city'), per]
    outdata.round(1).to_csv(outpath, mode='a', header=False)

    # mid (linwood)
    with open(outpath, 'a') as f:
        f.write('Mid aquifer nitrate-N (mg/L)\n')
    outdata = pd.DataFrame(index=scenarios, columns=percentiles)
    for scen, per in itertools.product(scenarios, percentiles):
        outdata.loc[scen, per] = data.loc[('Linwood', 'full_city'), (scen, per)] + chch_8.loc[
            ('Linwood', 'full_city'), per]
    outdata.round(1).to_csv(outpath, mode='a', header=False)

    # deep ( deep unnamed)
    with open(outpath, 'a') as f:
        f.write('Deep aquifer nitrate-N (mg/L)\n')
    outdata = pd.DataFrame(index=scenarios, columns=percentiles)
    for scen, per in itertools.product(scenarios, percentiles):
        outdata.loc[scen, per] = data.loc[('deep unnamed', 'full_city'), (scen, per)] + chch_8.loc[
            ('deep unnamed', 'full_city'), per]
    outdata.round(1).to_csv(outpath, mode='a', header=False)


def output_reduction_percentages(data_path, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pc5_mults = get_pc5pa_mults()
    data = pd.read_csv(data_path, index_col=0, header=[0, 1])
    current_measured = gen_well_targets('current_measured')
    current_measured.update(gen_stream_targets('current_measured'))

    # gmp
    gmp = data['gmp'].loc[:, ['5%', '50%', '95%']]

    # current pathways
    cpath = (gmp.transpose() * (1 + pc5_mults.loc[gmp.index] / 2)).transpose()

    # streams
    cols = ['GMP_50%', 'GMP_95%', 'cpath_50%', 'cpath_95%']
    for stream in streams:
        stream_long_targets = [6.9, current_measured[stream], 3.8, 1.0, 2.4]
        outdata = pd.DataFrame(index=[str(e) for e in stream_long_targets],
                               columns=cols, dtype=str)
        for tar in stream_long_targets:
            outdata.loc[str(tar), 'GMP_50%'] = '{}%'.format(int(round((1 - tar / gmp.loc[stream, '50%']) * 100, 0)))
            outdata.loc[str(tar), 'GMP_95%'] = '{}%'.format(int(round((1 - tar / gmp.loc[stream, '95%']) * 100, 0)))
            outdata.loc[str(tar), 'cpath_50%'] = '{}%'.format(int(round((1 - tar / cpath.loc[stream, '50%']) * 100, 0)))
            outdata.loc[str(tar), 'cpath_95%'] = '{}%'.format(int(round((1 - tar / cpath.loc[stream, '95%']) * 100, 0)))
        outdata.to_csv(os.path.join(outdir, '{}_reductions.csv'.format(stream)))

    # wdc_wells
    wdc_long_targets = [[5.6 for e in wdc_wells],
                        [8.5 for e in private_wells],
                        [10 for e in private_wells],
                        [current_measured[e] for e in wdc_wells]]
    for targets, name in zip(wdc_long_targets, ['half_mav',
                                                '3_4 mav',
                                                '90_mav',
                                                'current_measured']):
        outdata = pd.DataFrame(index=wdc_wells, columns=cols, dtype=str)
        for well, tar in zip(wdc_wells, targets):
            outdata.loc[well, 'GMP_50%'] = '{}%'.format(int(round((1 - tar / gmp.loc[well, '50%']) * 100, 0)))
            outdata.loc[well, 'GMP_95%'] = '{}%'.format(int(round((1 - tar / gmp.loc[well, '95%']) * 100, 0)))
            outdata.loc[well, 'cpath_50%'] = '{}%'.format(int(round((1 - tar / cpath.loc[well, '50%']) * 100, 0)))
            outdata.loc[well, 'cpath_95%'] = '{}%'.format(int(round((1 - tar / cpath.loc[well, '95%']) * 100, 0)))
        outdata.to_csv(os.path.join(outdir, 'wdc_wells_{}.csv'.format(name)))

    # private wells
    private_long_targets = [
        [5.6 for e in private_wells],
        [7.1 for e in private_wells],
        [8.5 for e in private_wells],
        [current_measured[e] for e in private_wells]]
    for targets, name in zip(private_long_targets, ['half_mav',
                                                    'no_exceed',
                                                    '3_4_mav',
                                                    'current_measured']):
        outdata = pd.DataFrame(index=private_wells, columns=cols, dtype=str)
        for well, tar in zip(private_wells, targets):
            outdata.loc[well, 'GMP_50%'] = '{}%'.format(int(round((1 - tar / gmp.loc[well, '50%']) * 100, 0)))
            outdata.loc[well, 'GMP_95%'] = '{}%'.format(int(round((1 - tar / gmp.loc[well, '95%']) * 100, 0)))
            outdata.loc[well, 'cpath_50%'] = '{}%'.format(int(round((1 - tar / cpath.loc[well, '50%']) * 100, 0)))
            outdata.loc[well, 'cpath_95%'] = '{}%'.format(int(round((1 - tar / cpath.loc[well, '95%']) * 100, 0)))
        outdata.round().to_csv(os.path.join(outdir, 'private_wells_{}.csv'.format(name)))


def output_area_weighted_reductions(well_path, wdc_wells_shps, outpath):
    """

    :param well_path:
    :param wdc_wells_shps: boolean true if wdc wells
    :param outpath:
    :return:
    """

    if wdc_wells_shps:
        private_shapes_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\wdc_use_mix"
    else:
        private_shapes_dir = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\private_wells_90_named_right"
    private_data = pd.read_csv(well_path, index_col=0)
    outdata = {}
    for key in private_data.keys():
        temp_arrays = []
        for well in private_data.index:
            if wdc_wells_shps:
                use_well = well.replace('wdc_', '')
            else:
                use_well = well
            temp = smt.shape_file_to_model_array(os.path.join(private_shapes_dir, '{}.shp'.format(use_well)), 'Id',
                                                 alltouched=True)
            temp[np.isfinite(temp)] = float(private_data.loc[well, key].replace('%', ''))
            temp_arrays.append(temp[np.newaxis, :, :])
        temp_arrays = np.concatenate(temp_arrays, axis=0)
        outdata[key] = {'05%': '{}%'.format(int(round(np.nanpercentile(temp_arrays, 5), 0))),
                        '25%': '{}%'.format(int(round(np.nanpercentile(temp_arrays, 25), 0))),
                        '50%': '{}%'.format(int(round(np.nanpercentile(temp_arrays, 50), 0))),
                        '75%': '{}%'.format(int(round(np.nanpercentile(temp_arrays, 75), 0))),
                        '95%': '{}%'.format(int(round(np.nanpercentile(temp_arrays, 95), 0))),
                        'max': '{}%'.format(int(round(np.nanmax(temp_arrays), 0))),
                        }
    outdata = pd.DataFrame(outdata)
    outdata.to_csv(outpath)


if __name__ == '__main__':
    datapaths = [
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\long_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv"
    ]
    prefixes = ['long_red_', 'all_scens_']
    if False:
        print('making current pathways table')
        output_current_pathways_table(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info")  # todo

        for data_path, prefix in zip(datapaths, prefixes):
            print('making future scenarios table for {}'.format(prefix))
            output_future_scenarios_table_without('all', data_path,
                                                  outdir=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info",
                                                  prefix=prefix)
            output_future_scenarios_table_withpa('all', data_path,
                                                 outdir=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info",
                                                 prefix=prefix)

    if False:
        print('making interzone tables for all scen')
        output_interzone_table(scenarios='all',
                               data_path=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\interzone\all_n_interzone.csv",
                               outpath=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\interzone_summary_all_scen.csv")
        print('making interzone tables for long scen')
        output_interzone_table(scenarios='all',
                               data_path=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\long_scens\interzone\all_n_interzone.csv",
                               outpath=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\interzone_summary_long_red.csv")
    if False:
        print('making reduction percentages')
        output_reduction_percentages(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\reductions")

    if True:
        print('making percentage reduction summaries')
        paths = glob(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\reductions\*_wells_*.csv")
        outdir = os.path.join(os.path.dirname(paths[0]), 'summaries')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for path in paths:
            wdc = 'wdc' in os.path.basename(path)
            outpath = os.path.join(outdir, 'summary_' + os.path.basename(path))
            output_area_weighted_reductions(path, wdc, outpath)
