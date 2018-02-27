# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/11/2017 4:13 PM
"""

from __future__ import division
from core import env
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    data = pd.read_csv(r"T:\Temp\temp_gw_files\test_emma_phi.csv")
    for key in data.keys():
        fig, ax = plt.subplots(figsize=(18.5, 9.5))
        ax.hist(data[key].dropna(),bins=100)
        ax.set_title(key)
    plt.show()
