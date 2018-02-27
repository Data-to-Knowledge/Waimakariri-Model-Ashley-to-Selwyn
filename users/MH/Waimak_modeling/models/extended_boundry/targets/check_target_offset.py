# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/07/2017 12:58 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
def calc_target_offset():
    all_targets = pd.read_csv(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/first_pass_head_targets_2008.csv"),
                              index_col=0)

    all_targets = all_targets.loc[:, [u'nztmx', u'nztmy', u'depth', u'rl_from_dem',
                                      u'ground_level', u'mid_screen_depth', u'mid_screen_elv',
                                      u'h2o_dpth_mean', u'readings', u'count_dry', u'count_flowing',
                                      u'h2o_elv_mean', u'owner_measured', u'bad_data', u'monitoring_well', u'samp_time_var',
                                      u'total_error_m', u'include_non-gap',
                                      u'm_error', u'users_error', u'scorrection_error', u'dem_error',
                                      u'low_rd_error']
                  ]

    vert_targets = pd.read_excel(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/Vertical gradient targets updated.xlsx"),
                                 sheetname='data_for_python', index_col=0)

    all_targets.loc[:,'vert_group'] = np.nan

    elv = smt.calc_elv_db()
    for i in vert_targets.index:
        all_targets.loc[i,'vert_group'] = vert_targets.loc[i,'group']
    number_of_values = len(all_targets.index)
    for num,i in enumerate(all_targets.index):
        if num%100 == 0:
            print ('completed {} of {}'.format(num,number_of_values))
        try:
            layer, row, col = smt.convert_coords_to_matix(all_targets.loc[i,'nztmx'],all_targets.loc[i,'nztmy'],all_targets.loc[i,'mid_screen_elv'])
            all_targets.loc[i,'layer'] = layer
            all_targets.loc[i,'row'] = row
            all_targets.loc[i,'col'] = col
            mx, my, mz = smt.convert_matrix_to_coords(row,col,layer,elv)
            all_targets.loc[i,'mx'] = mx
            all_targets.loc[i,'my'] = my
            all_targets.loc[i,'mz'] = mz
        except AssertionError as val:
            print(val)

    all_targets.to_csv(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_offsets.csv"))

def plot_values(data):
    data['z_offset'] = data.loc[:,'mid_screen_elv'] - data.loc[:,'mz']

    fig,(ax1,ax2) = plt.subplots(2,1, figsize=(18.5, 9.5))
    s1 = ax1.scatter(x=data.depth,y=data.z_offset)
    ax1.set_ylabel('z offset')
    ax1.set_xlabel('depth')

    s2 = ax2.scatter(x=data.ground_level,y=data.z_offset)
    ax2.set_ylabel('z offset')
    ax2.set_xlabel('ground level')

    fig2, axs = plt.subplots(3,4, figsize=(18.5, 9.5))
    vmin,vmax = -20,20
    model_xs, model_ys = smt.get_model_x_y()
    no_flow = smt.get_no_flow()
    for i in range(11):
        tempdata = data.loc[data.layer==i]
        ax = axs.flatten()[i]
        ax.set_aspect('equal')
        ax.contour(model_xs, model_ys, no_flow[i])
        ax.scatter(x=tempdata.nztmx,y=tempdata.nztmy,s=tempdata.depth,c=tempdata.z_offset,vmin=vmin, vmax=vmax)

    ax = axs.flatten()[11]
    ax.set_aspect('equal')
    ax.contour(model_xs, model_ys, no_flow[0])
    s = ax.scatter(x=data.nztmx,y=data.nztmy,c=data.z_offset,vmin=vmin, vmax=vmax)
    cb = fig2.colorbar(s)
    cb.set_label('zoffset')

    return fig,fig2


if __name__ == '__main__':
    calc_target_offset() # already run and it is quite slow
    data = pd.read_csv(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_offsets.csv"))
    print(data[data['include_non-gap']].readings.min())
    fig, fig2 = plot_values(data[data['include_non-gap']])

    print('done')
