# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/12/2017 10:08 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import flopy

if __name__ == '__main__':
    mask = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\raw_sw_samp_points\sfr\other\ashley_garry.shp".format(smt.sdp),'k',True)
    mask = np.isfinite(mask)
    idxs = smt.model_where(mask)
    data = nc.Dataset(r"C:\mh_waimak_model_data\post_filter1_budget.nc")

    outdata = [data.variables['stream leakage'][:,e[0],e[1]][:,np.newaxis] for e in idxs]
    outdata = np.concatenate(outdata,axis=1).sum(axis=1)/86400
    outdata = pd.Series(outdata,index=data.variables['nsmc_num'][:])

    mask = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\raw_sw_samp_points\sfr\other\cust_benit_cust.shp".format(smt.sdp),'k',True)
    mask = np.isfinite(mask)
    idxs = smt.model_where(mask)

    outdata_cust = [data.variables['stream leakage'][:,e[0],e[1]][:,np.newaxis] for e in idxs]
    outdata_cust = np.concatenate(outdata_cust,axis=1).sum(axis=1)/86400
    outdata_cust = pd.Series(outdata_cust,index=data.variables['nsmc_num'][:])

    cons = nc.Dataset(r"C:\mh_waimak_model_data\mednload_ucn.nc")
    chch_n_average = smt.shape_file_to_model_array(r"{}\m_ex_bd_inputs\shp\rough_chch.shp".format(smt.sdp),'Id',True)
    mask = np.isfinite(chch_n_average)
    idxs = smt.model_where(mask)
    out_chch_n_average = np.concatenate([np.array(cons.variables['mednload'][:, :, e[0], e[1]])[:,:,np.newaxis] for e in idxs],axis=2)
    out_chch_n_average = np.nanmean(out_chch_n_average,axis=2)
    out_chch_n_average = np.nanmean(out_chch_n_average,axis=1)
    out_chch_n_average = pd.Series(out_chch_n_average,index=cons.variables['nsmc_num'[:]])

    out = pd.merge(outdata.to_frame('ash'), out_chch_n_average.to_frame('con'),right_index=True,left_index=True)
    out = pd.merge(out,outdata_cust.to_frame('cust'),right_index=True,left_index=True)
    out.to_csv(r"C:\Users\matth\Downloads\ashley_losses_and_chch_con.csv")


