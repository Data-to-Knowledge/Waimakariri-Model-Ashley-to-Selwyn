# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 10/01/2018 5:45 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from core.stats.LR_class import LR

if __name__ == '__main__':
    bio = pd.read_csv(r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\ISO.CRDS.40")
    geo = pd.read_csv(r"C:\Users\MattH\Desktop\stuff_for_geothermal_paper\ISO.CRDS.115")

    bio_lr = LR(bio.CO2INV, bio.Delta_Raw_iCO2)
    geo_lr = LR(geo.CO2INV, geo.Delta_Raw_iCO2)

    fig, ax = plt.subplots(1)
    fig.suptitle('Keeling Plots')

    for lr, nm, c, y in zip([bio_lr,geo_lr],['Biogenic Signature', 'Geothermal Signature'],['k','grey'],[-20,-2.5]):
        ax.scatter(lr.x,lr.y, label=nm, color=c)
        temp_x = [0,0.001]
        ax.plot(temp_x,lr.predict(temp_x), color=c)
        ax.set_ylabel(u'$\delta^{13}C-CO_{2}$ (‰)')
        ax.set_xlabel(u'$1/[CO^{2}]$ $(ppmv^{-1})$')
        ax.set_ylim(-28,0)
        ax.set_xlim(temp_x)
        ax.text(0.00015,y,u'intercept: {} ‰'.format(round(lr.predict(0),1)),color=c)
    ax.legend()

    plt.show()
