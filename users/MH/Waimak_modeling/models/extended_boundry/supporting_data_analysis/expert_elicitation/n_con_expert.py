# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 12/12/2017 9:13 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import os

if __name__ == '__main__':
    base_emma_dir = r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells_rerun"
    experts = pd.read_excel(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\InterzoneTransferExpertJudgementWorkshop\Processed elcitation results.xlsx",
        sheetname='Distribution parameters')
    temp = pd.read_excel(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\InterzoneTransferExpertJudgementWorkshop\Processed elcitation results.xlsx",
        sheetname='EMMA')
    group_ids = temp.loc[:,['Target','Group_']].set_index('Group_')
    # combine all the expert PDFs and convert to normal distribution????
    outdata = pd.DataFrame(index=list(set(experts.group)),
                           columns=['mean', 'sd', '1st', '5th', '25th', '50th', '75th', '95th', '99th'])
    for group in set(experts.group):

        print(group)
        temp = []
        for expert, alpha, beta in experts.loc[experts.group == group, ['expert', 'alpha', 'beta']].itertuples(False,
                                                                                                               None):
            temp.append(np.random.beta(alpha, beta, 10000))
        expert_data = np.concatenate(temp)

        n = 100000
        # load EMMA results
        temp = []
        for target in set(group_ids.loc[int(group[-1]), 'Target']):
            temp.append(np.loadtxt(os.path.join(base_emma_dir, '{}_river.txt'.format(target))))


        emma_results = 1 - np.concatenate(temp)
        # eyre dilution
        eyre = np.random.normal(0.75, 0.15, n)
        emma_temp = np.random.choice(emma_results, n)
        expert_temp = np.random.choice(expert_data, n)

        # multiply by N load
        final_dist = emma_temp * expert_temp * eyre * 16.3


        # output
        outdata.loc[group, 'mean'] = final_dist.mean()
        outdata.loc[group, 'sd'] = final_dist.std()
        outdata.loc[group, '1st'] = np.percentile(final_dist, 1)
        outdata.loc[group, '5th'] = np.percentile(final_dist, 5)
        outdata.loc[group, '25th'] = np.percentile(final_dist, 25)
        outdata.loc[group, '50th'] = np.percentile(final_dist, 50)
        outdata.loc[group, '75th'] = np.percentile(final_dist, 75)
        outdata.loc[group, '95th'] = np.percentile(final_dist, 95)
        outdata.loc[group, '99th'] = np.percentile(final_dist, 99)

    outdata.to_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\InterzoneTransferExpertJudgementWorkshop\monte_carlo_expert_n.csv")
