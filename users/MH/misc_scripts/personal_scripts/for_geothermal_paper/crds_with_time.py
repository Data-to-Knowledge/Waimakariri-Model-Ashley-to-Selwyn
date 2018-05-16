# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/05/2018 9:28 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from core.stats.LR_class import LR
from glob import glob
import os

# make crds plots with time... 2 subplots on ontop of the other one showing d13c and the other showing [co2]

if __name__ == '__main__':
    data_paths = glob(r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\ISO_files\ISO_files\ISO.CRDS.*")
    outdir = r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\all_crds_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data = []
    for path in data_paths:
        geo = pd.read_csv(path)
        data.append(geo.X12CO2.values)
    data = np.concatenate(data).flatten()

    for path in data_paths:
        site = os.path.basename(path).split('.')[-1]
        geo = pd.read_csv(path)
        geo.loc[:, 'time'] = geo.loc[:, 'EPOCH_TIME'] - geo.EPOCH_TIME.min()
        fig, (ax1, ax2) = plt.subplots(2)
        time_lims = [geo.time.min(), geo.time.max()]
        ax1.set_xlim(time_lims)
        ax2.set_xlim(time_lims)
        ax2.set_xlabel('Time (s)')

        ax1.set_ylabel(u'$[CO^{2}]$ $(ppmv^{-1})$')
        ax1.set_ylim([300, geo.X12CO2.max()])

        ax2.set_ylabel(u'$\delta^{13}C-CO_{2}$ (‰)')
        ax2.set_ylim(-28, 0)
        time = geo.time
        # plot co2 con ax1
        con = geo.X12CO2
        ax1.scatter(time, con, color='k')

        # plot d13c con ax2
        d13c = geo.Delta_Raw_iCO2
        ax2.scatter(time, d13c, color='k')

        fig.suptitle('chamber evolution for site {}'.format(site))
        fig.savefig(os.path.join(outdir, 'chamber_evol_site_{:03d}.png'.format(int(site))))
        plt.close(fig)

