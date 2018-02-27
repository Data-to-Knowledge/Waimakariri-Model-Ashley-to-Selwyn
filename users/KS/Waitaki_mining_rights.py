# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:49:27 2017

@author: KateSt
"""

#### Import hydro class
from core.classes.hydro import hydro
import numpy as np
import pandas as pd

mtypes1 = 'flow'
mtypes2 = 'flow_m'

sites1 = [71106, 71102, 71167, 71178, 71105, 71104]
sites2 = [1841, 2367, 1711435, 71105, 71111, 171159, 171158, 171147]
sites3 = [71106, 71178]

qual_codes = qual_codes=[50, 40, 30, 20, 10, 11, 21, 18]

recorder_data = hydro().get_data(mtypes='flow', sites=sites1, qual_codes=qual_codes) #adds data to h1

awakino_data = hydro().rd_csv(csv_path=r'S:\Surface Water\shared\projects\Waitaki Mining Rights\71105_awakino.csv', sites=None, mtypes='flow', values=None, dformat='wide', time='Time', multiindex=False, skiprows=0)


#### Get out gauging data for 71106 and 71178 combined record

recorder_gaugings = hydro().get_data(mtypes='flow_m', sites=sites1, qual_codes=qual_codes)

stats_df = recorder_gaugings._base_stats

stats_df2 = stats_df.reset_index()
stats_df3 = stats_df2.set_index('site')

start_time_71178 = stats_df3.loc[71178, 'start_time']
end_time_71178 = stats_df3.loc[71178, 'end_time']

start_time_71167 = stats_df3.loc[71167, 'start_time']
end_time_71167 = stats_df3.loc[71167, 'end_time']

df = recorder_gaugings.sel_ts(mtypes=mtypes2, pivot=False)
df = df.unstack(level=0)
df = df.unstack(level=0)

df.loc[df2.index < start_time_71178, ('flow_m', 71178)] = df.loc[df.index < start_time_71178, ('flow_m', 71167)]
df.loc[df2.index < start_time_71167, ('flow_m', 71178)] = df.loc[df.index < start_time_71167, ('flow_m', 71102)]

drop_old_sites = df.drop([('flow_m', 71102), ('flow_m', 71167)], axis=1)

reformat = drop_old_sites.stack(['mtype', 'site'])
reformat.name = 'values'

reformat_1 = reformat.reset_index()

#### get gauging data for other sites

gauging_site_data = hydro().get_data(mtypes=mtypes2, sites=sites2, qual_codes=qual_codes)

all_site_gaugings = gauging_site_data.add_data(data=reformat_1, time='time', values='values', mtypes='mtype', sites='site', dformat='long') #adds dataframe to hydroclass object

export_csv1 = r'S:\Surface Water\shared\projects\Waitaki Mining Rights\all_site_gaugings.csv'


####

all_site_data = recorder_data.get_data(mtypes=mtypes2, sites=sites2, qual_codes=qual_codes)

all_site_data_2 = all_site_data.add_data(data=reformat_1, time='time', values='values', mtypes='mtype', sites='site', dformat='long')

all_site_data_3 = all_site_data_2.sel_ts(pivot=True)

# Calculate median flow at primary sites
medians = all_site_data_3.median(axis=0)
median_71106 = medians.loc['flow',71106]
median_71178 = medians.loc['flow',71178]

# Calcualte upper quartile flow at primary sites
uqs = all_site_data_3.quantile(q=0.75, axis=0)
uq_71106 = uqs.loc['flow',71106]
uq_71178 = uqs.loc['flow',71178]

### replace flow with gaugings where
idx = (~all_site_data_3['flow_m'][71106].isnull())
all_site_data_3.loc[idx, ('flow', 71106)] = all_site_data_3.loc[idx, ('flow_m', 71106)]

idx1 = (~all_site_data_3['flow_m'][71178].isnull())
all_site_data_3.loc[idx1, ('flow', 71178)] = all_site_data_3.loc[idx1, ('flow_m', 71178)]
###


idx2 = np.logical_or(all_site_data_3['flow'][71106] <= uq_71106, all_site_data_3['flow'][71178] <= uq_71178)

all_site_data_4 = all_site_data_3.loc[idx2]

all_site_data_5 = all_site_data_4.drop([('flow_m', 71106), ('flow_m', 71178)], axis=1)

all_site_data_drop = all_site_data_5['flow_m'].dropna(axis=0, how='all')

all_site_data_6 = all_site_data_5.loc[all_site_data_drop.index]

export_csv1 = r'S:\Surface Water\temp\Waitaki_Mining_Rights_1.csv'

all_site_data_6.to_csv(export_csv1)

####
#drop_recorder_gaugings = all_site_data_5.drop([('flow_m', 71106), ('flow_m', 71178)], axis=1)

reformat2 = all_site_data.stack(['mtype', 'site'])
reformat2.name = 'values'

reformat3 = reformat2.reset_index()
####

all_site_data_6 = hydro().add_data(data=reformat3, time='time', values='values', mtypes='mtype', sites='site', dformat='long')

all_site_data_6.sel_ts(pivot=True)

export_csv1 = r'S:\Surface Water\temp\Waitaki_Mining_Rights_1.csv'

all_site_data_6.to_csv(export_csv1, pivot=True)

####

new1, reg = all_site_data_6.flow_reg(y=sites2, x=sites3, p_val=0.5, min_yrs=5, min_obs=5)
malf6 = new1.malf7d()
new1
reg
malf6



####
"""
all_site_data_3['flow_all'][71106] = all_site_data_3['flow_m'][71106]
all_site_data_3['flow_all'][71178] = all_site_data_3['flow_m'][71178]

idx = (all_site_data_3['flow_m'][71106].isnull())
all_site_data_3.loc[idx, ('flow_all', 71106)] = all_site_data_3.loc[idx, ('flow', 71106)]

idx1 = (all_site_data_3['flow_m'][71178].isnull())
all_site_data_3.loc[idx1, ('flow_all', 71178)] = all_site_data_3.loc[idx, ('flow', 71178)]

all_site_data_4 = all_site_data_3.drop([('flow', 71106), ('flow', 71178)], axis=1)

all_site_data_6 = all_site_data_3['flow_m'].dropna(axis=0, how='all')

all_site_data_6 = all_site_data_4.loc[all_site_data_5.index]
"""
