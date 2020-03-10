# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 4/05/2018 8:54 AM
"""

from __future__ import division
import env
import pandas as pd


def get_pc5pa_mults(inc_interzone=False):
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

    if inc_interzone:
        interzone_pas = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulatio"
                                            r"ns and results\ex_bd_va\n_results\interzone\pa_rules_interzone\load_over"
                                            r"views_with_paN.csv"), index_col=0)
        interzone_pas = interzone_pas.loc[:, 'total_pa_N_kg'] / interzone_pas.loc[:, 'gmp_nload_kg']
        interzone_pas.name = 'pa_mult'
        outdata = pd.concat((outdata,interzone_pas))

    return outdata


def get_mt3d_current_pathway_n(add_pc5=True, inc_interzone=False):
    """

    :param add_pc5: boolean if True add the additional concentration associated with the pc5 pa rules
    :param inc_interzone: boolean if True add the interzone targets
    :return:
    """
    pa_mults = get_pc5pa_mults(inc_interzone=inc_interzone)
    gmp_results = pd.read_csv(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations an"
                                      r"d results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\a"
                                      r"ll_n_waimak_zone.csv"), index_col=0, header=[0, 1])['gmp']

    if inc_interzone:
        interzone_results = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulation"
                                        r"s and results\ex_bd_va\zc_n_sols\all_scens\interzone\all_n_interzone.csv",
                                        index_col=[0,1], header=[0,1])['gmp']
        gmp_results.loc['conservative_interzone'] = interzone_results.loc[('deep unnamed', 'full_city')]
        gmp_results.loc['highly_likely_interzone'] = interzone_results.loc[('deep unnamed', 'full_city')]

    test = (1 + pa_mults.loc[gmp_results.index]) #todo should this pamult be divided by 2?
    if not add_pc5:
        test.loc[:] = 1

    out = (gmp_results.transpose() * test).transpose()
    return out



if __name__ == '__main__':
    test = get_pc5pa_mults()
    test = get_mt3d_current_pathway_n(inc_interzone=True)
    print('done')
