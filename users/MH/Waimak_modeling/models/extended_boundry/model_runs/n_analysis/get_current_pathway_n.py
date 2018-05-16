# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/05/2018 8:54 AM
"""

from __future__ import division
from core import env
import pandas as pd


def get_pc5pa_mults():
    wdc_pas = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations an"
                                  r"d results\ex_bd_va\n_results\wdc_use_mix\pa_rules_wdc_use_mix\load_overviews_wi"
                                  r"th_paN.csv"), index_col=0)
    wdc_pas = wdc_pas.loc[:, 'total_pa_N_kg'] / wdc_pas.loc[:, 'gmp_nload_kg']
    wdc_pas.index = ['wdc_{}'.format(e) for e in wdc_pas.index]
    wdc_pas.name = 'pa_mult'

    private_pas = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations a"
                                      r"nd results\ex_bd_va\n_results\private_wells_90\pa_rules_private_wells_90\load_o"
                                      r"verviews_with_paN.csv"), index_col=0)
    private_pas = private_pas.loc[:, 'total_pa_N_kg'] / private_pas.loc[:, 'gmp_nload_kg']
    private_pas.name = 'pa_mult'

    stream_pas = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations an"
                                     r"d results\ex_bd_va\n_results\nwaimak_springfeds\pa_rules_nwaimak_springfeds\loa"
                                     r"d_overviews_with_paN.csv"), index_col=0)
    stream_pas = stream_pas.loc[:, 'total_pa_N_kg'] / stream_pas.loc[:, 'gmp_nload_kg']
    stream_pas.name = 'pa_mult'

    outdata = pd.concat((wdc_pas, private_pas, stream_pas))

    return outdata


def get_mt3d_current_pathway_n():
    pa_mults = get_pc5pa_mults()
    gmp_results = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations an"
                                      r"d results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\a"
                                      r"ll_n_waimak_zone.csv"), index_col=0, header=[0, 1])['gmp']
    test = (1 + pa_mults.loc[gmp_results.index])
    out = (gmp_results.transpose() * test).transpose()
    return out



if __name__ == '__main__':
    test = get_pc5pa_mults()
    print('done')
