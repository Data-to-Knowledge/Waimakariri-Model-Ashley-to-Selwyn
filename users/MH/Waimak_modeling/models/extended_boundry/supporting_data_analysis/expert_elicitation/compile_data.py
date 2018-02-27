# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/08/2017 9:42 AM
"""

from __future__ import division
from core import env
import pandas as pd
import os
from glob import glob
from file_set_up import user_codes, users, minor_qoi, major_qoi, values, minor_qoi_shapes
from warnings import warn
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm, norm

all_dir = "C:/Users/MattH/OneDrive - Environment Canterbury/waimakariri_elicitation/all_qois"
if not os.path.exists(all_dir):
    os.makedirs(all_dir)
def combine_QOI (qoi):
    qoi_dir = '{}/{}'.format(all_dir,qoi)
    if not os.path.exists(qoi_dir):
        os.makedirs(qoi_dir)
    allout_data = []
    for sub_qoi in minor_qoi[qoi]:
        out_data = pd.DataFrame(index=users.sort(),columns=values)
        for user in users:
            path = "C:\Users\MattH\OneDrive - Environment Canterbury\waimakariri_elicitation\{u}\W_elicitation_21-08-2017_{u}_{q}.xlsx".format(u=user,q=qoi)
            temp = pd.read_excel(path, index_col=0)
            out_data.loc[user,values] = temp.loc[sub_qoi,values]

        if out_data.isnull().any().any():
            raise ValueError('{} has null values from users: {}'.format(sub_qoi, list(out_data.index[out_data.isnull().any(1)])))

        out_data.to_csv('{}/{}_combined.csv'.format(qoi_dir,sub_qoi))
        allout_data.append(out_data)
    return allout_data

def plot_pdfs(datasets, names, colors, dims, linear_pool=False, distribution=skewnorm):
    if not isinstance(datasets, tuple):
        raise ValueError('expected tuple of datasets')
    fig, axs = plt.subplots(*dims)
    out_axes = {}
    for i, data in enumerate(datasets):
        ax = np.atleast_1d(axs).flatten()[i]
        x = np.linspace(data.lower_limit.min(), data.upper_limit.max(), 100)
        for j, person in enumerate(data.index):
            _min, _q1, m, _q2, _max = data.loc[person, ['lower_limit', 'lower_quartile', 'median',
                                                        'upper_quartile','upper_limit']]

            samples = 10000

            q1 = np.linspace(_min, _q1, samples / 4)
            q2 = np.linspace(_q1, m, samples / 4)
            q3 = np.linspace(m, _q2, samples / 4)
            q4 = np.linspace(_q2, _max, samples / 4)

            temp = np.concatenate((q1, q2, q3, q4))
            if linear_pool:
                if j==0:
                    alldata=temp
                else:
                    alldata = np.concatenate((alldata,temp))
                line_style = '--'
            else:
                line_style = '-'

            params = distribution.fit(temp)
            ax.plot(x, distribution.pdf(x, *params), linestyle=line_style, color=colors[j], lw=1, label=person)
        if linear_pool:
            params = distribution.fit(alldata)
            x = np.linspace(distribution.ppf(0.001, *params),
                            distribution.ppf(0.999, *params))
            ax.plot(x, distribution.pdf(x, *params), linestyle='-', color='blue', lw=2,
                    label='linear pool')
            print ('parameters for {}: {}'.format(names[i], params))
        ax.set_title(names[i])
        out_axes[names[i]] = ax
    plt.legend()
    return fig, out_axes


def combine_and_plot(qoi,linear_pool=False,additional=False, distribution=skewnorm):
    colors = ['red', 'darkgreen', 'fuchsia', 'black','orange','teal','darkmagenta']
    datasets = tuple(combine_QOI(qoi))
    if additional:
        linear_pool = True
    fig, axs = plot_pdfs(datasets,minor_qoi[qoi],colors,minor_qoi_shapes[qoi],linear_pool,distribution=distribution)

    if additional:
        for name in minor_qoi[qoi]:
            ax = axs[name]
            data = pd.read_excel(r"C:\Users\MattH\OneDrive - Environment Canterbury\waimakariri_elicitation\group\W_elicitation_21-08-2017_group_{}.xlsx".format(qoi))
            x = np.linspace(data.loc[name,'lower_limit'], data.loc[name, 'upper_limit'], 100)
            _min, _q1, m, _q2, _max = data.loc[name, ['lower_limit', 'lower_quartile', 'median',
                                                            'upper_quartile', 'upper_limit']]
            samples = 10000

            q1 = np.linspace(_min, _q1, samples / 4)
            q2 = np.linspace(_q1, m, samples / 4)
            q3 = np.linspace(m, _q2, samples / 4)
            q4 = np.linspace(_q2, _max, samples / 4)

            temp = np.concatenate((q1, q2, q3, q4))
            line_style = '-'

            params = distribution.fit(temp)
            ax.plot(x, distribution.pdf(x, *params), linestyle=line_style, color='green', lw=3,
                    label='group')

        plt.legend()

    plt.show(fig)


if __name__ == '__main__':
    print 'start'
    qois = ['LSR', 'pumping', 'LRZ_flux', 'offshore', 'possible_ks', 'race_loss']
    combine_and_plot(major_qoi[0],linear_pool=True,additional=False,distribution=norm)






