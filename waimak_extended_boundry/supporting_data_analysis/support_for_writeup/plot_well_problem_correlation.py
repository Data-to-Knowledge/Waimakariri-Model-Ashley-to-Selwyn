# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/10/2018 9:18 AM
"""

from __future__ import division
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_runs.n_analysis.interzone_n import get_chch_area_zones
import pandas as pd
import numpy as np
import os

zone_labels = {'central_chch':(1.56961e6, 5.18287e6),
               'north_western_tla':(1.55654e6,5.18821e6),
               'western_chch':(1.56123e6,5.17747e6),
               'north_east':(1.57029e6,5.19122e6),
               'southern_chch':(1.56978e6,5.17679e6)}

if __name__ == '__main__':
    outpath = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\model_checks\pumpin_problem"
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    data = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\model_checks\pumpin_problem\correlations.csv",index_col=0)
    chch_zones = get_chch_area_zones()
    plot_zones = ['central_chch', 'north_western_tla', 'western_chch', 'north_east', 'southern_chch']
    zone_numbers = smt.get_empty_model_grid() * np.nan
    for i,zone in enumerate(plot_zones):
        zone_numbers[chch_zones[zone]] = i+1

    for layer in range(10):
        fig, ax = smt.plt_matrix(zone_numbers,base_map=True,cmap='Dark2')
        ax.set_xlim([1.551E6, 1.58204e6])
        ax.set_ylim([5.17239E6, 5.19694e6])


        for lab, coord in zone_labels.items():
            key = 'layer_{:02d}_zone_{}'.format(layer, lab)
            text = '{}\n cor: {:.3e}\n p: {:.3e}'.format(lab.replace('_',' '),
                                                 data.loc[key,'pearsonr'],
                                                 data.loc[key,'pearsonp'])
            ax.text(coord[0],coord[1],text,bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})

        fig.savefig(os.path.join(outpath,'mapping_layer_{:02d}.png'.format(layer)))

    print('done')
