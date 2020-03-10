"""
Author: matth
Date Created: 27/03/2018 9:09 AM
"""

from __future__ import division
import shutil
from glob import glob
import os

if __name__ == '__main__':
    paths = glob("D:\mh_waimak_models\stocastic_forward\*_results")
    dests = [os.path.join(r"K:\mh_modeling\stocastic_forward", os.path.basename(e)) for e in paths]
    for path, dest in zip(paths, dests):
        shutil.copytree(path, dest)
