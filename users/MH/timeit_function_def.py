from __future__ import division
from users.MH.Waimak_modeling.models.extended_boundry.m_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.cwms_index import get_zone_array_index
from scipy.interpolate import rbf
import netCDF4 as nc
import numpy as np
from core import env
import shutil
import os
import pandas as pd
import flopy
import itertools
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.extract_data import open_path_file_as_df
base_data = np.random.rand(50000000)
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.n_analysis.interzone_n import _np_describe

def timeit_function():
    _np_describe(base_data)
def timeit_function2():
    _np_describe(np.random.choice(base_data,100000))



if __name__ == '__main__':
    print 'hello world'
    bot = smt.calc_elv_db()[1]
    hds = flopy.utils.HeadFile(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\supporting_data_for_scripts\ex_bd_va_sdp\from_gns\NsmcBase\AW20171024_2_i2_optver\i2\mf_aw_ex.hds").get_alldata()[0][0]
    idx = get_zone_array_index('waimak')
    hds[hds>1e20] = np.nan
    dif1 = hds-bot
    dif1[dif1<0] = np.nan
    hds[~idx] = np.nan
    dif2 = hds-bot
    dif2[dif2<0] = np.nan

    for nm, data in zip(['all','waimak'],[dif1,dif2]):
        print nm
        print 'min: {}'.format(np.nanmin(data))
        print '1 :{}'.format(np.nanpercentile(data,1))
        print '5 :{}'.format(np.nanpercentile(data,5))
        print '25 :{}'.format(np.nanpercentile(data,25))
        print '50 :{}'.format(np.nanpercentile(data,50))
        print '75 :{}'.format(np.nanpercentile(data,75))
        print '95 :{}'.format(np.nanpercentile(data,95))
        print '99 :{}'.format(np.nanpercentile(data,99))
        print 'max: {}'.format(np.nanmax(data))
        print 'u :{}'.format(np.nanmean(data))
        print 'std: {}'.format(np.nanstd(data))



    print 'done'