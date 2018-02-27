# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/08/2017 1:46 PM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt, _get_constant_heads
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.colors import from_levels_and_colors, LogNorm
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.targets.gen_target_arrays import get_head_targets
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.drn_packages import _get_drn_spd




outdir = "{}/figs_for_report".format(smt.sdp)

carpets = ['d_ash_c','d_chch_c','d_dlin_c','d_dsel_c','d_ulin_c', 'd_usel_c','d_cust_c']
drn_data = _get_drn_spd(1,1)
temp = drn_data.loc[~np.in1d(drn_data.parameter_group,carpets)]
temp = smt.add_mxmy_to_df(temp)
temp.to_csv('{}/drn_data_non_carpet.csv'.format(outdir))

maper = dict(zip(carpets,range(1,len(carpets)+1)))
drn_data = drn_data.loc[np.in1d(drn_data.parameter_group,carpets)]
drn_data = drn_data.replace({'parameter_group': maper})
smt.array_to_raster('{}/carpet_drains2.tif'.format(outdir),smt.df_to_array(drn_data,'parameter_group'))




cmap, norm = from_levels_and_colors([-1, 0, 1, 2], ['blue', 'black', 'white'])

# plots of head targets
heads = get_head_targets()
for layer in range(10):
    print layer
    temp = heads.loc[heads.layer==layer]
    no_flow = smt.get_no_flow(layer)
    fig, ax = smt.plt_matrix(no_flow,cmap=cmap,color_bar=False, title='head targets for layer {}'.format(layer+1),
                             no_flow_layer=None,norm=norm)
    sc = ax.scatter(temp.nztmx,temp.nztmy,c=temp.weight,cmap='plasma', norm=LogNorm())
    fig.colorbar(sc, extend='max')
    fig.set_size_inches(5.5, 5.5)
    fig.savefig('{}/head_targets_layer{:02d}'.format(outdir,layer+1))

# te wai and constant heads plot
cheads = _get_constant_heads()[0]
te_wai = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/te_waihora.shp".format(smt.sdp), 'ID', True)
cheads[np.isfinite(te_wai)] = 1.5
fig,ax = smt.plt_matrix(cheads)
fig.set_size_inches(5.5, 5.5)
fig.savefig('{}/cheads.png'.format(outdir))


khpaths = glob(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\InitialParamaters\pilot_pts_2-8-2017\KH_ppk_*.txt"))

# pilot points plot
for path in khpaths:
    data = pd.read_table(path,names=['site','x','y','layer','val'])
    layer = int(path.split('_')[-1].split('.')[0]) -1
    no_flow = smt.get_no_flow(layer)
    fig,ax = smt.plt_matrix(no_flow,title='KV and KH pilot points layer {}'.format(layer+1),no_flow_layer=None,cmap=cmap,
                            color_bar=False,norm=norm)
    ax.scatter(data.x,data.y,c='red')
    fig.set_size_inches(5.5, 5.5)
    fig.savefig('{}/pilot_points_layer{:02d}.png'.format(outdir,layer))

