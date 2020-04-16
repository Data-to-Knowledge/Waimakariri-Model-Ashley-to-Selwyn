# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 25/10/2017 3:09 PM
"""

from __future__ import division
import os

#todo look through documentation

def trans_readme(data_dir,out_dir):
    inpath = os.path.join(data_dir,'READ_ME.txt')
    outpath = os.path.join(out_dir,'READ_ME.txt')
    if os.path.exists(inpath):
        with open(inpath) as f:
            lines = f.readlines()
        with open(outpath,'w') as f:
            f.writelines(lines)

