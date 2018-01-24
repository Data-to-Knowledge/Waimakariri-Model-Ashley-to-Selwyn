from __future__ import division
import numpy as np
import pandas as pd
import flopy
import glob
import matplotlib.pyplot as plt
import os
from core.ecan_io import rd_sql, sql_db
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import get_wel_spd
from pykrige.ok import OrdinaryKriging as okrig
import geopandas as gpd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.extract_data import open_path_file_as_df

missing_wells = ['M35/0043',
                 'M35/0575',
                 'M35/11576',
                 'M35/11798',
                 'M35/11865',
                 'M35/11866',
                 'M35/5247',
                 'M35/6456',
                 'M35/6473',
                 'M35/6538',
                 'M35/7254',
                 'M35/9283',
                 'M35/9360',
                 'M35/9362',
                 'M35/9363',
                 'M35/9375']
path_path = r"D:\mh_waimak_models\modpath_reverse_base\test_multiple\weak\NsmcReal000005\NsmcReal000005_weak.mppth"
group_mapper_path = r"D:\mh_waimak_models\modpath_reverse_base\test_multiple\weak\NsmcReal000005\NsmcReal000005_weak_group_mapper.csv"
group_mapper = pd.read_csv(group_mapper_path, index_col=0, names=['key', 'val'])['val'].to_dict()

data = open_path_file_as_df(path_path)
data.rename(columns={'Layer':'layer','Row':'row', 'Column':'col'}, inplace=True)
data.replace({'Particle_Group':group_mapper},inplace=True)
data = data.loc[np.in1d(data.Particle_Group, missing_wells)]
print 'adding mx,my'
data = smt.add_mxmy_to_df(data)
print 'saving data'
data.to_csv(r"T:\Temp\temp_gw_files\missing_particles.csv")
print('done')


