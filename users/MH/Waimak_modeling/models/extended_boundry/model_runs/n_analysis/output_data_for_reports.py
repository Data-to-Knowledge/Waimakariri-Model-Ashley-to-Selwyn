# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/05/2018 8:57 AM
"""

from __future__ import division
from core import env
import os
import pandas as pd
import numpy as np
import itertools
from get_current_pathway_n import get_pc5pa_mults
from percentage_reduction_maps import wdc_wells, private_wells, streams, gen_stream_targets, gen_well_targets


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
        temp = ['{} ({}-{})'.format(round(m,2), round(f,2), round(n,2)) for
                f, m, n in data['cmp'].loc[receptor_set, ['5%', '50%', '95%']].itertuples(False, None)]
        outdata.loc[:, 'cmp'] = temp

        # gmp
        gmp = data['gmp'].loc[receptor_set, ['5%', '50%', '95%']]
        temp = ['{} ({}-{})'.format(round(m,2), round(f,2), round(n,2)) for f, m, n in gmp.itertuples(False, None)]
        outdata.loc[:, 'gmp'] = temp

        # pc5pa
        pc5 = (gmp.transpose() * (1 + pc5_mults.loc[gmp.index])).transpose()
        temp = ['{} ({}-{})'.format(round(m,2), round(f,2), round(n,2)) for f, m, n in pc5.itertuples(False, None)]
        outdata.loc[:, 'pc5pa'] = temp

        # cpath
        cpath = (gmp.transpose() * (1 + pc5_mults.loc[gmp.index] / 2)).transpose()
        temp = ['{} ({}-{})'.format(round(m,2), round(f,2), round(n,2)) for f, m, n in cpath.itertuples(False, None)]
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
            temp2 = ['{} ({}-{})'.format(round(m, 2), round(f, 2), round(n, 2)) for f, m, n in
                    temp.itertuples(False, None)]
            outdata.loc[:, scen] = temp2
        outdata.to_csv(outpath)



def output_future_scenarios_table_without(scenarios, data_path, outdir,prefix=''):
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
            temp2 = ['{} ({}-{})'.format(round(m, 2), round(f, 2), round(n, 2)) for f, m, n in
                     temp.itertuples(False, None)]
            outdata.loc[:, scen] = temp2
        outdata.to_csv(outpath)

def output_interzone_table(scenarios, data_path, outpath):
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
    outdata.to_csv(outpath, mode='a', header=False)

    # mid (linwood)
    with open(outpath, 'a') as f:
        f.write('Mid aquifer nitrate-N (mg/L)\n')
    outdata = pd.DataFrame(index=scenarios, columns=percentiles)
    for scen, per in itertools.product(scenarios, percentiles):
        outdata.loc[scen, per] = data.loc[('Linwood', 'full_city'), (scen, per)] + chch_8.loc[
            ('Linwood', 'full_city'), per]
    outdata.to_csv(outpath, mode='a', header=False)

    # deep ( deep unnamed)
    with open(outpath, 'a') as f:
        f.write('Deep aquifer nitrate-N (mg/L)\n')
    outdata = pd.DataFrame(index=scenarios, columns=percentiles)
    for scen, per in itertools.product(scenarios, percentiles):
        outdata.loc[scen, per] = data.loc[('deep unnamed', 'full_city'), (scen, per)] + chch_8.loc[
            ('deep unnamed', 'full_city'), per]
    outdata.to_csv(outpath, mode='a', header=False)


if __name__ == '__main__':
    datapaths = [
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\long_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
        r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv"
    ]
    prefixes = ['long_red_', 'all_scens_']
    if True:
        output_current_pathways_table(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\all_n_waimak_zone.csv",
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info") #todo

        for data_path,prefix in zip(datapaths,prefixes):
            output_future_scenarios_table_without('all',data_path,
                                                  outdir=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info",
                                                  prefix=prefix)
            output_future_scenarios_table_withpa('all',data_path,
                                                 outdir=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info",
                                                 prefix=prefix)

    if True:
        output_interzone_table(scenarios='all',
                               data_path=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\interzone\all_n_interzone.csv",
                               outpath=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\interzone_summary_all_scen.csv")
        output_interzone_table(scenarios='all',
                               data_path=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\long_scens\interzone\all_n_interzone.csv",
                               outpath=r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\summary_info\interzone_summary_long_red.csv")
