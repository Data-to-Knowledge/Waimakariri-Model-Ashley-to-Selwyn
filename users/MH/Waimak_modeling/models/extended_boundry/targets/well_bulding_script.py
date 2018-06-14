# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 25/07/2017 11:13 AM
"""

import pandas as pd
import numpy as np
import flopy_mh as flopy

#this script was passed to brioch for inclusion in the pest optimisation process
# the influx wells are presently set to 1 m3/s so the muliplier can range between 0 and 5
well_data = pd.read_csv()# set path
muliplier = 1 # set muliplier
well_data.loc[well_data.type=='lr_boundry_flux','flux'] *= muliplier
g = well_data.groupby(['layer', 'row', 'col'])
outdata = g.aggregate({'flux': np.sum}).reset_index()
outdata = outdata.rename(columns={'layer': 'k', 'row': 'i', 'col': 'j'}).to_records(False)
outdata = outdata.astype(flopy.modflow.ModflowWel.get_default_dtype())

#write into file