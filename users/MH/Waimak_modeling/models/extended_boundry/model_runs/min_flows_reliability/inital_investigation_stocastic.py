# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/03/2018 2:02 PM
"""

from __future__ import division
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    [u'', u'', u'', u'cam_marshes_s',
     u'', u'kaiapoi_harpers_s', u'',
     u'', u'', u'saltwater_factory_s',
     u'southbrook_marshes_s', u'taranaki_gressons_s', u'taranaki_preeces_s',
     u'waikuku_sh1_s', u'ash_ash_est', u'ash_est_all', u'cam_end_s',
     u'', u'taranaki_end_s', u'',
     u'kaiapoi_end']

    mapper = {u'cam_youngs': 'cam_bramleys_s',
              u'courtenay_neeves': 'courtenay_kaiapoi_s',
              u'northbrook_marsh': 'northbrook_marshes_s',
              u'ohoka_island': 'ohoka_island_s',
              u'saltwater_toppings': 'saltwater_end_s',
              u'silverstream_neeves': 'kaiapoi_island_s',
              u'southbrook_marsh': 'southbrook_marshes_s',
              u'taranaki_preeces': 'taranaki_preeces_s',
              u'waikuku_waikuku-beach-rd': 'waikuku_end_s',
              u'ashley_sh1': 'ashley_sh1',
              u'cust_oxford': 'cust_skewbridge',
              u'cust_threlkelds': 'cust_skewbridge'}
    for key in mapper.keys():
        flow_data = pd.read_csv(
            r"K:\mh_modeling\stocastic_forward\condensed_data\raw_data\{}_raw_data.csv".format(key), index_col=0)
        wai_data = pd.read_csv(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_str_data.csv",
            index_col=0)

        data = pd.merge(flow_data, pd.DataFrame(wai_data.loc[:, mapper[key]]), right_index=True, left_index=True)
        data = data.rename(columns={mapper[key]: 'wai_per'})
        fig, ax = plt.subplots(figsize=(18.5, 9.5))

        colours = ['black',
                   'red',
                   'gold',
                   'darkgreen',
                   'navy',
                   'darkorchid',
                   'm',
                   'darkorange',
                   'maroon',
                   'blue',
                   'green',
                   'grey',
                   'lawngreen',
                   'cyan',
                   'purple',
                   'deepskyblue', ] #todo get more useful colors (need 15)
        for i, d in enumerate(set(data.keys()) - {'wai_per'}):
            ax.scatter(data.loc[:, 'wai_per'], data.loc[:, d], label=d, c=colours[i])
            ax.set_xlabel('wai_per')
            ax.set_ylabel('flow change')
            ax.set_title(key)
            ax.legend()

    plt.show()
    print('done')
    # looking through this work it appears that there is no real corellation between waimak percentage and flow changes.
    # just use as standard