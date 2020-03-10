"""
Author: matth
Date Created: 24/05/2017 9:34 AM
"""

from __future__ import division
import netCDF4 as nc
import pandas as pd
import numpy as np
from users.MH.Waimak_modeling.supporting_data_path import sdp


sites_org = pd.read_csv("{}/inputs/vs_nums.txt".format(sdp))

sites = sites_org['Network']

vcsn = nc.Dataset(r"S:\Data\VirtualClimate\vcsn_precip_et_2016-06-06.nc")

v_sites = np.array(nc.chartostring(vcsn.variables['site'][:]))
time = nc.num2date(vcsn.variables['time'][:],vcsn.variables['time'].units)
et = vcsn.variables['ET'][:][:,np.in1d(v_sites,sites)]
precip = vcsn.variables['precip'][:][:,np.in1d(v_sites,sites)]

precip_at_sites = np.nanmean(precip, axis=0)
precip_at_sites_st = np.nanstd(precip, axis=0)

et_average = np.nanmean(et, axis=1)
precip_average = precip.mean(axis=1)
out_data = pd.DataFrame({'time': time, 'et': et_average, 'precip': precip_average})
out_data.to_csv("{}/inputs/spatially_mean_vscn.csv".format(sdp))

all_et = pd.DataFrame(data=et,columns=sites)

all_precip = pd.DataFrame(data=precip,columns=sites)
all_et.to_csv("{}/inputs/all_et_vscn.csv".format(sdp))
all_precip.to_csv("{}/inputs/all_precip_vscn.csv".format(sdp))


