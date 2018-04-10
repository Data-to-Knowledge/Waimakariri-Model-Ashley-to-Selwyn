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

def timeit_function():
    with open(r"C:\Users\MattH\Downloads\test_mt3d\test") as f:
        num_lines = sum([1 for l in f])
    data = pd.read_table(r"C:\Users\MattH\Downloads\test_mt3d\test", skiprows=int(num_lines - 2*(num_lines-2)/1228),
                         delim_whitespace=True,
                         names=[u'TIME', u'SFR-NODE', u'SFR-CONCENTRATION', u'FLOWGW', u'GW-CONC'],
                         dtype={u'TIME':float, u'SFR-NODE':int, u'SFR-CONCENTRATION':float, u'FLOWGW':float, u'GW-CONC':float},
                         low_memory=True)
    data = data.loc[np.isclose(data.TIME, data.TIME.max())]
def timeit_function2():
    data = pd.read_table(r"C:\Users\MattH\Downloads\test_mt3d\test", skiprows=1, delim_whitespace=True,
                         low_memory=True)
    data = data.loc[np.isclose(data.TIME, data.TIME.max())]



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