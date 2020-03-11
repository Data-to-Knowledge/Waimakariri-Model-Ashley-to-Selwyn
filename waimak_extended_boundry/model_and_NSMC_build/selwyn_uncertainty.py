# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 14/08/2017 12:57 PM
"""

from __future__ import division
from core.classes.hydro import hydro
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


h1 = hydro().get_data(['flow','flow_m'], sites=[167,68001],from_date='2008-01-01')
h2 = hydro().get_data(['flow','flow_m'], sites=[167,68001])

h1.data = h1.data.loc[h1.data<10]

plt.hist(h1.data.loc['flow',68001],bins=100)
plt.axvline(h1.data.loc['flow',68001].mean(),color='r')
plt.axvline(h1.data.loc['flow',68001].median(),color='b')

plt.close()

white = pd.DataFrame(h2.data.loc['flow',68001])
coal = pd.DataFrame(h2.data.loc['flow_m',167])
temp = pd.merge(white,coal,left_index=True,right_index=True)

regr = stats.linregress(temp.data_x.values,temp.data_y.values)

# Make predictions using the testing set
y_pred = temp.data_x* regr[0] + regr[1]


plt.scatter(temp.data_x,temp.data_y)
plt.plot(temp.data_x,y_pred)
plt.text(0.5,8,'m = {:.2f}x + {:.2f}'.format(regr[0],regr[1]))
plt.ylabel('coalgate')
plt.xlabel('whitecliffs')
plt.title('Selwyn River Flows')

print 'done'