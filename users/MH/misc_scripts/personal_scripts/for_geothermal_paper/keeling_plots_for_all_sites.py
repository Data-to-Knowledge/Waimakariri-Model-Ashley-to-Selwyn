# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 2/05/2018 9:21 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from core.stats.LR_class import LR
from glob import glob
import os

if __name__ == '__main__':
    data_paths = glob(r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\ISO_files\ISO_files\ISO.CRDS.*")
    outdir=r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\all_keeling_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data = []
    for path in data_paths:
        geo = pd.read_csv(path)
        data.append(geo.CO2INV.values)
    data = np.concatenate(data).flatten()

    for path in data_paths:
        geo = pd.read_csv(path)
        lr = LR(geo.CO2INV, geo.Delta_Raw_iCO2)
        site = os.path.basename(path).split('.')[-1]
        fig, ax = plt.subplots()
        temp_x = [0,0.003]
        ax.scatter(lr.x, lr.y, color='k')
        ax.plot(temp_x,lr.predict(temp_x), color='k', linestyle='--')
        ax.set_ylabel(u'$\delta^{13}C-CO_{2}$ (‰)')
        ax.set_xlabel(u'$1/[CO^{2}]$ $(ppmv^{-1})$')
        ax.set_ylim(-28,0)
        ax.set_xlim(temp_x)
        ax.text(0.00015,-2,u'intercept: {} ‰'.format(round(lr.predict(0),1)))
        ax.set_title('keeling plot for site {}'.format(site))
        fig.savefig(os.path.join(outdir,'keeling_plot_site_{:03d}.png'.format(int(site))))
        plt.close(fig)
