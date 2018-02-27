# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 23/08/2017 12:49 PM
"""

from __future__ import division
from core import env
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skewnorm, norm


def individual_pdf(path, distribution):
    print('')
    print(path)
    data = pd.read_excel(path, index_col=0)

    for person in data.index:
        fig, ax = plt.subplots(1,1)
        _min, _q1, m, _q2, _max = data.loc[person, ['lower_limit', 'lower_quartile', 'median',
                                                    'upper_quartile', 'upper_limit']]

        samples = 10000

        q1 = np.linspace(_min, _q1, samples / 4)
        q2 = np.linspace(_q1, m, samples / 4)
        q3 = np.linspace(m, _q2, samples / 4)
        q4 = np.linspace(_q2, _max, samples / 4)

        temp = np.concatenate((q1, q2, q3, q4))

        params = distribution.fit(temp)
        print '{}: {}; percent error at 95th: {}'.format(person,params,params[1]/params[0]*200)
        x = np.linspace(_min, _max, 100)
        ax.plot(x, distribution.pdf(x, *params), label=person)
        #plt.show()
        plt.close(fig)


if __name__ == '__main__':
    paths = [r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Expert judgement elicitation\elicitation_23-08-2017_zeb_etheridge_ChchSFS.xlsx",
             r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Expert judgement elicitation\elicitation_23-08-2017_zeb_etheridge_AshleyRiverlosses.xlsx",
             r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Expert judgement elicitation\elicitation_22-08-2017_zeb_etheridge_racelosses.xlsx",
             r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Expert judgement elicitation\elicitation_22-08-2017_zeb_etheridge_northernboundary.xlsx"]
    for path in paths:
        individual_pdf(path, norm)
