"""
Author: matth
Date Created: 16/06/2017 2:37 PM
"""

from __future__ import division
import users.MH.Waimak_modeling.model_tools as mt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main(segments = None, elv_db = None):
    if elv_db is None:
        elv = mt.calc_elv_db()
    str_data = pd.DataFrame(mt.get_base_str())
    for val in str_data.index:
        k,i,j = str_data.loc[val,['k','i','j']]
        str_data.loc[val,'glevel'] = elv[k,i,j] # get the top elevation
    if segments is None:
        segments = set(str_data['segment'])
    else:
        segments = np.atleast_1d(segments)
    for seg in segments:
        temp = str_data[str_data['segment']==seg]
        start, stop = temp['stop'].iloc[0], temp['stop'].iloc[-1]
        fig, ax = plt.subplots(1,1, figsize=(18.5, 9.5))
        ax.plot(temp['reach'],temp['sbot'],label='bot',color='r')
        ax.plot(temp['reach'],temp['stop'],label='top',color='b')
        ax.plot(temp['reach'],temp['glevel'],label='ground level',color='k')
        ax.set_title(seg)
        ax.legend()
        plt.show()

if __name__ == '__main__':
    main()


