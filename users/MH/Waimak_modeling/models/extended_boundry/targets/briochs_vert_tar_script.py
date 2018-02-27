# -*- coding: utf-8 -*-
"""
Script to generated vertical targets from vertical_gradient_target.csv

Output obs2obs equations and pst observation data

Created on Mon Aug 21 11:07:43 2017

@author: briochh
"""

import os
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.targets.gen_target_arrays import get_vertical_gradient_targets

def get_vert_targets_full():
    vertwells = get_vertical_gradient_targets().reset_index()  # read in Matt H csv
    vertwells.rename(columns={'Well No': 'Well'}, inplace=True)  # change name with space in
    vertwells['Well'] = vertwells['Well'].str.replace('/',
                                                      '_')  # replace '/' in name with '_' - consitent with bash methods for wells
    # vertwells['weight']=vertwells['weight']/(vertwells['weight']**0.5) # if we want to convert to 1/sqrt(sigma) - actually do that later
    gb = vertwells.sort_values('layer').groupby(
        'group')  # sort first so that subsequent groups are in layer order and we can use indexes to effiiecntly calculate head diffs

    # new df groupnumber,upperlayer in diff, lowerlayer in diff,codename for diff (obsname),lower well name in diff, upper well name in diff, (lower-upper)diff, 'composite weight'
    vgroups = pd.DataFrame(columns=['group', 'ulay', 'llay', 'name', 'well', 'minus', 'hdif', 'cwgt'])
    # %%
    out_data = pd.DataFrame()
    c = 0  # index counter for new df
    # individual groups
    for group, g in gb:  # loop over groups
        # print(name)
        for i in range(g.shape[0] - 1):  # within each group
            #    print([g.layer.iloc[i],g.layer.iloc[i+1]])
            ulay = int(g.layer.iloc[i] + 1)  # upper layer
            llay = int(g.layer.iloc[i + 1] + 1)  # lower layer
            code = '{}_{}_{}'.format(g.code.iloc[i], llay, ulay)  # obs code name
            well = g.Well.iloc[i + 1]  # lower well
            minus = g.Well.iloc[i]  # upper well (subtracted from lower)
            # aggregated error (assuming well error is 1/sigma)
            # add sigmas and calsulate 1/sqrt(sum(sigmas))
            compwgt = 1 / (((1 / g.weight.iloc[i + 1])**2 + (1 / g.weight.iloc[i])**2) ** 0.5)
            # save to df
            vgroups.loc[c, ['group', 'ulay', 'llay', 'name', 'well', 'minus', 'hdif', 'cwgt']] = \
                [group, ulay, llay, code, well, minus, g.obs.iloc[i + 1] - g.obs.iloc[i], compwgt]
            c += 1  # increment to next entry in df
            # there are more than 2 entries for a group and we are currently looking at the bottom entry
            if g.shape[0] > 2 and i + 1 == g.shape[0] - 1:
                # head diff across top and bottom entries
                ulay = int(g.layer.iloc[0] + 1)  # top layer
                llay = int(g.layer.iloc[i + 1] + 1)  # bottom layer
                code = '{}_{}_{}'.format(g.code.iloc[i], llay, ulay)
                well = g.Well.iloc[i + 1]  # bottom well
                minus = g.Well.iloc[0]  # top well
                compwgt = 1 / (((1 / g.weight.iloc[i + 1])** 2 + (1 / g.weight.iloc[0])** 2) ** 0.5)
                vgroups.loc[c, ['group', 'ulay', 'llay', 'name', 'well', 'minus', 'hdif', 'cwgt']] = \
                    [group, ulay, llay, code, well, minus, g.obs.iloc[i + 1] - g.obs.iloc[0], compwgt]
                c += 1
                # print(np.diff(g.obs))
    # %% save results
    # equations file for obs2obs
    # instrcutions file for obs2obs
    # observations file for pst
            eqfmt = '{:11} = {:10} - {:10}\n'.format  # obsname = lowerwell - upper well
            insfmt = 'l1 [{}]22:50\n'.format  # l1 [obsname]22:49
            obsfmt = '{:{name_width}}{:10.6f}   {:10.6f}  vert\n'.format  # obsname  headdifference weight obsgroup
        out_data = pd.concat((out_data,vgroups)) # matt added
    return vgroups
if __name__ == '__main__':
    vgroups = get_vert_targets_full()

    print 'test'

