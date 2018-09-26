# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/10/2017 11:29 AM
"""

from __future__ import division
from core import env
import os
import pandas as pd

# a scipt to test whether there are any overlaps... looking for duplicates
def find_str_in_scripts(directory, string):
    all_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                all_paths.append(os.path.join(root,file))
    scripts = []
    pickle_paths = []
    for path in all_paths:
        script = os.path.basename(path)
        with open(path) as f:
            lines = f.readlines()
        for line in lines:
            if string in line:
                scripts.append(script)
                pickle_paths.append(line)
    data = pd.DataFrame({'lines':pickle_paths,'scripts':scripts})
    data.to_csv(os.path.join(env.temp("temp_gw_files"),'find_string_in_scripts.csv'))



if __name__ == '__main__':
    find_str_in_scripts(r'C:\Users\MattH\PycharmProjects\Ecan.Science.Python.Base',
                        'import')
