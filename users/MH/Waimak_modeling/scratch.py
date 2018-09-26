from __future__ import division
from glob import glob
import shutil
from core import env
import os
import pandas as pd

if __name__ == '__main__':

    strings = ' '.join(list(pd.read_csv(os.path.join(env.temp("temp_gw_files"),'find_string_in_scripts.csv'))['lines']))
