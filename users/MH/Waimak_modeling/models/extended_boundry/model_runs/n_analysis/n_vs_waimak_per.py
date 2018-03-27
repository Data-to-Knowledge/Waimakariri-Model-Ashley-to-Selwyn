# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/03/2018 3:16 PM
"""

from __future__ import division
from core import env
import pandas as pd
import matplotlib.pyplot as plt
from core.stats.LR_class import LR
from matplotlib.offsetbox import AnchoredText


if __name__ == '__main__':
    # saltwater increases with increasing alpine river water because the Northern boundary flux was mistakenly not set
    # in the LSR models.

    npath =r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points\raw_stocastic_set_str_data.csv"
    waipath =r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\waimak_per_results_at_points\raw_stocastic_set_str_data.csv"
    n = pd.read_csv(npath,index_col=0)
    sites = n.keys()
    wai = pd.read_csv(waipath, index_col=0)
    data = pd.merge(n,wai,left_index=True,right_index=True,suffixes=('_n','_wai'))
    for site in sites:
        if data.loc[:,'{}_wai'.format(site)].isnull().all():
            continue
        fig, ax = plt.subplots(figsize=(18.5, 9.5))
        model = LR(data.loc[:,'{}_wai'.format(site)],data.loc[:,'{}_n'.format(site)])
        ax.scatter(data.loc[:,'{}_wai'.format(site)],data.loc[:,'{}_n'.format(site)])
        ax.plot(data.loc[:,'{}_wai'.format(site)], model.predict(data.loc[:,'{}_wai'.format(site)]))
        anchored_text = AnchoredText("formula: {} \n R2: {}".format(model.formula, model.adj_rval), loc=2)
        ax.add_artist(anchored_text)
        ax.set_ylabel('N')
        ax.set_xlabel('alpine_river_fraction')
        ax.set_title(site)

    plt.show()
