from __future__ import division
import numpy as np
import pandas as pd
import geopandas as gpd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import matplotlib.pyplot as plt
from copy import deepcopy
from glob import glob
from core.classes.hydro import hydro, all_mtypes
from matplotlib.colors import from_levels_and_colors

khpaths = glob(r"C:\Users\MattH\Downloads\ppt_forplproc\*_ppk_*.txt")
cmap, norm = from_levels_and_colors([-1, 0, 1, 2], ['blue', 'black', 'white'])

for path in khpaths:
    data =pd.read_table(path,names=['site','x','y','layer','val'])
    layer = data.layer.iloc[0]-1
    no_flow = smt.get_no_flow(layer)
    fig,ax = smt.plt_matrix(no_flow,title='KV and KH pilot points layer {}'.format(layer),no_flow_layer=None,cmap=cmap,
                            norm=norm)
    ax.scatter(data.x,data.y)
    plt.show(fig)
