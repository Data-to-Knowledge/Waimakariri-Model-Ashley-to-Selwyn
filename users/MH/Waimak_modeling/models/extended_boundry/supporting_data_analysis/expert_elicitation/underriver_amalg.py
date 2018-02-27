# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/11/2017 11:38 AM
"""

from __future__ import division
from core import env
import pandas as pd
from glob import glob
import os

def get_match_data(path):
    with open(path) as f:
        lines = f.readlines()
    max = [float(e.split(' ')[-1].strip()) for e in lines if 'Minimum' in e][0]
    min = [float(e.split(' ')[-1].strip()) for e in lines if 'Maximum' in e][0]
    alpha = [float(e.split(' = ')[1].split(';')[0]) for e in lines if 'Scaled Beta:' in e][0]
    beta = [float(e.split(' = ')[2].split(';')[0]) for e in lines if 'Scaled Beta:' in e][0]
    norm_mean = [float(e.split(' = ')[1].split(';')[0]) for e in lines if 'Normal:' in e][0]
    norm_sd = [float(e.split(' = ')[2].split(';')[0]) for e in lines if 'Normal:' in e][0]

    return min,max,alpha,beta,norm_mean,norm_sd

def amalg_data(outpath):
    paths = glob(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quantity\ExpertJudgement_Files\All\*_*.txt")

    outdata = {'group':[],
               'expert': [],
               'min':[],
               'max':[],
               'alpha':[],
               'beta':[],
               'norm_mean':[],
               'norm_sd':[]}

    for path in paths:
        outdata['group'].append(os.path.basename(path).split('_')[0])
        outdata['expert'].append(os.path.basename(path).split('_')[1].lower().replace('.txt',''))

        min,max,alpha,beta,norm_mean,norm_sd = get_match_data(path)

        for key in set(outdata.keys())-{'group','expert'}:
            outdata[key].append(eval(key))

    outdata = pd.DataFrame(outdata)
    outdata.to_csv(outpath)

if __name__ == '__main__':
    amalg_data(r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quantity\ExpertJudgement_Files\All\amalgamated.csv")