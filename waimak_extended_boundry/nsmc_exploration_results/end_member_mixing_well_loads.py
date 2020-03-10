# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/11/2017 9:03 AM
"""

from __future__ import division
from waimak_extended_boundry import get_wel_spd
import numpy as np


def make_end_mixing_load_layers(outpath, version):
    """
    makes a pandas data frame for load layers for MT3d for the endmember mixing stuffs
    :param version: either 'inland' or 'alpine_river'
    :return:
    """

    wells = get_wel_spd(3)

    if version == 'inland':
        outwells = wells.loc[np.in1d(wells.type, ['river', 'boundry_flux'])]
    elif version == 'alpine_river':
        outwells = wells.loc[wells.type == 'race']
    else:
        raise ValueError('unexpected version entry {}'.format(version))
    outwells = outwells.loc[:, ['layer', 'row', 'col']].astype(int)
    outwells += 1
    outwells.loc[:, 'conc'] = 1.0
    outwells.loc[:, 'type'] = 2
    outwells.index.name='well'
    outwells.to_csv(outpath, sep=' ')

def make_median_n_load_layers(outpath):
    """
    makes a pandas data frame for load layers for MT3d for the endmember mixing stuffs
    :param version: either 'inland' or 'alpine_river'
    :return:
    """

    wells = get_wel_spd(3)

    outwells = wells.loc[np.in1d(wells.type, ['boundry_flux', 'race']) & (wells.zone=='n_wai')]
    outwells.loc[:, 'conc'] = 0.1
    outwells.loc[(outwells.type=='boundry_flux') & (outwells.col>150), 'conc'] = 2.0
    outwells = outwells.loc[:, ['layer', 'row', 'col','conc']]
    outwells.loc[:,['layer', 'row', 'col']] = outwells.loc[:,['layer', 'row', 'col']].astype(int)+ 1
    outwells.loc[:, 'type'] = 2
    outwells.index.name='well'
    outwells.to_csv(outpath, sep=' ')


if __name__ == '__main__':
    make_median_n_load_layers(r"C:\Users\Public\Desktop\median_n_load_wells.txt")
    make_end_mixing_load_layers(r"C:\Users\Public\Desktop\alpine_river_wells.txt", 'alpine_river')
    make_end_mixing_load_layers(r"C:\Users\Public\Desktop\inland_wells.txt", 'inland')
