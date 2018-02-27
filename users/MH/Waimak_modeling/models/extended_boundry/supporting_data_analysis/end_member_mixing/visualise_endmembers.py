# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/11/2017 8:50 AM
"""

from __future__ import division
from core import env
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    tracers = [
        'cl',
        'o18',
        'k'
    ]
    endmembers = [
        'inland',
        'coastal',
        'river',
        'eyre'
    ]

    means = {
        'cl': {
            'inland': 10.38,
            'coastal': 25.43,
            'river': 1.05,
            'eyre': 3.67
        },
        'o18': {
            'inland': -8.76,
            'coastal': -8.00,
            'river': -9.25,
            'eyre': -8.90
        },
        'k': {
            'inland': 1.23,
            'coastal': 1.67,
            'river': 0.57,
            'eyre': 1.17
        }
    }

    sds = {
        'cl': {
            'inland': 4.21,
            'coastal': 9.49,
            'river': 0.15,
            'eyre': 0.29
        },
        'o18': {
            'inland': 0.21,
            'coastal': 0.64,
            'river': 0.22,
            'eyre': 0.54
        },
        'k': {
            'inland': 0.31,
            'coastal': 0.30,
            'river': 0.15,
            'eyre': 0.12
        }
    }

    a_errors = {
        'cl': 0.5,
        'o18': 0.2,
        'k': 0.15
    }

    outdir =r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Groundwater Quality\End member mixing model\Additional target wells\4_members_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for trace in tracers:
        fig, ax = plt.subplots(figsize=(18.5, 9.5))
        for em in endmembers:
            data = np.random.normal(means[trace][em], sds[trace][em] + a_errors[trace], 100000)
            ax.hist(data, bins=101, label=em, alpha=0.5)

        fig.suptitle('tracer: {}'.format(trace))
        ax.legend()
        fig.savefig(os.path.join(outdir, '{}_endmembers.png'.format(trace)))


