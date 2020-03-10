# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/12/2017 3:28 PM
"""

from __future__ import division
import flopy_mh as flopy
import pandas as pd
from waimak_extended_boundry import smt
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # this was a quick test... the endpoints are similar for the weak source instance, but for strong source, it's terrible
    strong = pd.DataFrame(flopy.utils.EndpointFile(r"C:\Users\MattH\Desktop\test_reverse_modpath_strong\test_reverse.mpend").get_alldata())
    strong = strong.loc[strong['k'] == 0, ['i', 'j']].values
    weak = pd.DataFrame(flopy.utils.EndpointFile(r"C:\Users\MattH\Desktop\test_reverse_modpath_weak\test_reverse.mpend").get_alldata())
    weak = weak.loc[weak['k'] == 0, ['i', 'j']].values
    strong_out = smt.get_empty_model_grid().astype(bool)
    strong_out[strong[:,0],strong[:,1]] = True
    weak_out = smt.get_empty_model_grid().astype(bool)
    weak_out[weak[:,0],weak[:,1]] = True
    smt.plt_matrix(weak_out,title='weak_endpoint')
    smt.plt_matrix(strong_out,title='strong_endpoint')
    plt.show()



