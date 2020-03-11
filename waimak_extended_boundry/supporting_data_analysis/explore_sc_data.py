# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 21/07/2017 10:41 AM
depreciated due to changes in core
"""

from __future__ import division
import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from waimak_extended_boundry import smt
from core.ecan_io import rd_sql, sql_db

if __name__ == '__main__':
    all_wells = pd.read_csv('{}/all_wells_row_col_layer.csv'.format(smt.sdp),index_col=0)

    data =rd_sql(col_names=['WELL_NO','SpecificCapacity'], **sql_db.wells_db.well_details)
    data.loc[:,'WELL_NO'] = data.WELL_NO.str.replace('_','/')
    data = data.set_index('WELL_NO')
    data = data.rename(columns={'SPECIFIC_CAPACITY':'sc'})
    data = data.loc[(data.SpecificCapacity.notnull()) & (data.SpecificCapacity>0)]

    data = pd.merge(data,all_wells, how='inner', left_index=True, right_index=True)

    data.loc[:,'lnsc'] = np.log(data.SpecificCapacity)
    no_flow = smt.get_no_flow().astype(bool).astype(int)
    model_x, model_y = smt.get_model_x_y()



    for i in range(smt.layers):
        temp = data.loc[data.layer==i]
        fig,ax = plt.subplots(1,1,figsize=(18.5, 9.5))
        ax.contour(model_x,model_y, no_flow[i])
        s = ax.scatter(temp.nztmx,temp.nztmy,c=temp.lnsc, vmin=data.lnsc.min(),vmax=data.lnsc.max())
        fig.colorbar(s)
        ax.set_title('layer {}'.format(i))
        ax.set_aspect('equal')
        ax.set_xlim(model_x.min(),model_x.max())
        ax.set_ylim(model_y.min(),model_y.max())
        fig.savefig(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Model Grid\model_layering\sc_data_plots\layer{:02d}".format(i)))



    print 'done'