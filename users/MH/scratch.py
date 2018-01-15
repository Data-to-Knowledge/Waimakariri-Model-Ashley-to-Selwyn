from __future__ import division
from users.MH.Waimak_modeling.models.extended_boundry.m_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import create_wel_package
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from scipy.interpolate import rbf
import netCDF4 as nc
import numpy as np
from core import env
import shutil
import os
import itertools
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.extract_data import open_path_file_as_df

path_path = r"C:\Users\MattH\Desktop\test_reverse_modpath_strong\test_reverse.mppth"

drop_names = [
    'Time_Point_Index',
    'Cumulative_Time_Step',
    'Tracking_Time',
    'Global_X',
    'Global_Y',
    'Global_Z',
    'Grid',
    'Local_X',
    'Local_Y',
    'Local_Z',
    'Line_Segment_Index',
]
print('reading data')

data = open_path_file_as_df(path_path)
print('simplifying data')
data.drop(drop_names, 1, inplace=True)

sites = list(itertools.product([1],range(50,150),range(50,160)))

def timeit_function():
    temp_data = data.set_index(['Layer','Row','Column']).loc[sites].reset_index()
    return temp_data
def timeit_function2():
    temp_data = ['{}-{}-{}'.format(l,r,c) for l,r,c in data[['Layer','Row','Column']].itertuples(False,None)]
    sites2 = ['{}-{}-{}'.format(l,r,c) for l,r,c in sites]
    idx = np.in1d(temp_data, sites2)
    out_data = data.loc[idx]
    return out_data




if __name__ == '__main__':
    print 'hello world'
    temp = timeit_function()
    temp2 = timeit_function2()

    print 'done'