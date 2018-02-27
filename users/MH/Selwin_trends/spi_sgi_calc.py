"""
Author: matth
Date Created: 20/04/2017 8:49 AM
"""

from __future__ import division
from core import env
from core.stats.standard_precip.spi import SPI
import pandas as pd
import numpy as np
#this is not done and we need to.
input_data = pd.read_table(r"T:\Temp\Fouad\SPI_Calculator_input_output\input.txt", delim_whitespace=True)
spi = SPI()
spi.set_rolling_window_params(center=False)
spi.set_distribution_params('gamma')
data = spi.calculate(np.array(input_data['rain']))

print (data[0:10])
