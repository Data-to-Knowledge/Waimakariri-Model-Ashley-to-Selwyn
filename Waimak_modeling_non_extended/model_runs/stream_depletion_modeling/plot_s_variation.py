"""
Author: matth
Date Created: 1/06/2017 2:11 PM
"""

from __future__ import division

import matplotlib.pyplot as plt

import env
from s_cal_targets import get_s_cal_data
from users.MH.Waimak_modeling.model_tools.get_str_flows import create_hydrographs_for_transient

if __name__ == '__main__':
    stream_data = {}
    stream_data['real'] = get_s_cal_data()

    runs = ["D:\mattH\python_wm_runs\s_explore_transient\perc_25_ss_1.63783413356e-05_sy_0.0401643660801",
            r"D:\mattH\python_wm_runs\s_explore_transient\perc_50_ss_0.000169836616135_sy_0.000472232863272",
            r"D:\mattH\python_wm_runs\s_explore_transient\perc_75_ss_5.29569944766e-05_sy_0.0111526462969",
            r"D:\mattH\python_wm_runs\s_explore_transient\perc_90_ss_2.74998610469e-05_sy_0.0328226289808"]
    percentiles = [25, 50, 75, 90]

    for per, run in zip(percentiles,runs):
        stream_data[per] = create_hydrographs_for_transient(run)

    fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize=(18.5, 9.5))
    for site, ax in zip(['kaiapoi_harpers', 'cam_youngs', 'cust_threlkelds'], [ax1,ax2,ax3]):
        ax.set_title(site)
        ax.set_ylabel('m3/day')
        ax.set_xlabel('month')
        for run, c in zip(['real',25, 50,75,90],['k','r','y','b','purple']):
            x = stream_data[run].index
            y = stream_data[run][site]
            if isinstance(run,str):
                lab = run
            else:
                lab = 'S: {} percentile'.format(run)
            ax.plot(x, y, color=c, label=lab)
            ax.legend(loc=1)

    fig.tight_layout()
    fig.savefig(env.sci(r'Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\River and stream boundaries and targets\Transient data'))
