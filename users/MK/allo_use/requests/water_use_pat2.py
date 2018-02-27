# -*- coding: utf-8 -*-
"""
Examples on how to use HydroPandas.
"""

#### Import hydro class
from core.classes.hydro import hydro, all_mtypes

#################################################
#### Load data

### Parameters

mtypes6 = 'usage'
sites4 = ['K37/0419', 'K37/0494', 'K37/2170', 'K37/2171', 'K37/2391', 'K37/2502', 'K37/2706', 'K37/2707', 'K37/3205', 'K37/3206', 'K37/3207', 'K38/1512', 'K38/1513']

### From the MSSQL server (the easy way) - Loads in both the time series data and the geo locations

use1 = hydro().get_data(mtypes=mtypes6, sites=sites4)


#############################################
#### Saving data

export_csv1 = r'E:\ecan\local\Projects\requests\patrick\2017-07-19\water_use.csv'

### Export a time series csv
use2 = use1.sel_ts(pivot=True)
use3 = use2.round(1)

use3.to_csv(export_csv1)






