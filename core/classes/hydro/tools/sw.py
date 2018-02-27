# -*- coding: utf-8 -*-
"""
Surface water tools for the hydro class.
"""

from pandas import DataFrame, concat, Series
from core.ts.sw.stats import malf7d as malf_fun
from core.ts.sw.stats import malf_reg as malf_reg_fun
from core.ts.sw.stats import flow_reg as flow_reg_fun

#######################################
#### Adapted functions

def malf7d(self, sites=None, w_month='JUN', max_missing=90, malf_min=0.9, intervals=[10, 20, 30, 40], return_alfs=False, num_years=False, export_path=None, export_name_malf='malf.csv', export_name_alf='alf.csv', export_name_mis='alf_missing_data.csv'):

    data = self.sel_ts(mtypes='river_flow_cont', sites=sites, pivot=True)
    if data.index.inferred_freq != 'D':
        data = data.resample('D').mean()
    malf_set = malf_fun(data, w_month, max_missing, malf_min, intervals, return_alfs, num_years, export_path, export_name_malf, export_name_alf, export_name_mis)
    return(malf_set)



#def malf_reg(self, y, x=None, buffer_dis=None, min_yrs=10, min_obs=10, w_month='JUN', max_missing=120, malf_min=0.9, intervals=[10, 20, 30, 40], return_alfs=False):
#
#    #### Initial check
#    if not isinstance(y, list):
#        raise ValueError('y must be a list of site numbers')
#
#    #### Remove data based on restriction parameters
#
#    ### recorder flow removal
#    if 'flow' in self.mtypes:
#        fstats = self.stats(mtypes='flow')
#        sites = fstats[fstats['Tot data yrs'] >= min_yrs].index.tolist()
#    else:
#        raise ValueError('Load some flow data in!')
#
#    ### gauging flow removal
#    if 'm_flow' in self.mtypes:
#        m_count = self.sel_ts('m_flow', y).groupby(level='site').count()
#        m_sites = m_count[m_count >= min_obs].index.tolist()
##        sites.extend(m_sites)
#
#    #### Create new hydro object with updated sites
#    hp1 = self.sel(mtypes='flow', sites=sites)
#    hp2 = self.sel(mtypes='m_flow', sites=m_sites)
#
#    hydro1 = hp1.combine(hp2)
#
#    #### Extract/create comparison dictionary
#    if x is None:
#        if buffer_dis is None:
#            raise ValueError('If x is None, then buffer_dis should be defined.')
#        hydro2 = hydro1._comp_by_buffer(buffer_dis)
#        comp_dict = {i: list(hydro2.comp_dict[i]) for i in m_sites}
#    elif isinstance(x, list):
#        comp_dict = {i: x for i in m_sites}
#
#    #### Run regression
#
##    malf_lst = []
##    for j in comp_dict:
##        x_data = hydro2.sel_ts(mtypes='flow', sites=comp_dict[j], pivot=True)['flow']
##        y_data = hydro2.sel_ts(mtypes='m_flow', sites=j).loc(axis=0)['m_flow', j]
##        y_data.name = j
##        y_data = DataFrame(y_data)
##        reg1, malf1 = malf_reg_fun(x_data, y_data, roll=False, min_yrs=min_yrs, min_obs=min_obs, w_month=w_month, max_missing=max_missing, malf_min=malf_min, intervals=intervals)[0]
##        malf_lst.append(malf1)
##        reg_lst.append(reg1)
##    malf = concat(malf_lst)
#
#    malf = DataFrame()
#    reg = DataFrame()
#    for j in comp_dict:
#        x_data = hydro2.sel_ts(mtypes='flow', sites=comp_dict[j], pivot=True)
#        y_data = hydro2.sel_ts(mtypes='m_flow', sites=j, pivot=True, resample=False)
#        y_data.name = j
#        y_data = DataFrame(y_data)
#        reg1, y_new1 = flow_reg(x_data, y_data, min_obs=min_obs, p_val=0.05, logs=False, make_ts=True)
#        malf1 = malf_fun(y_new1, w_month, max_missing, malf_min, intervals, return_alfs)
#        malf = concat([malf, malf1])
#        reg = concat([reg, reg1])
#    return(malf, reg)


def flow_reg(self, y, x=None, y_mtype='river_flow_disc_qc', x_mtype='river_flow_cont_qc', buffer_dis=None, min_yrs=10, min_obs=10, logs=False, p_val=0.05, below_median=False):
    """
    Function to do a simple linear regression between gauging sites and recorder sites.
    The output is a hydro class object and the regression info.
    """
    #### Initial check
    if not isinstance(y, list):
        raise ValueError('y must be a list of site numbers')

    #### Remove data based on restriction parameters

    if x_mtype == 'river_flow_cont_qc':
        ### recorder flow removal
        if 'river_flow_cont_qc' in self.mtypes:
            fstats = self.stats(mtypes='river_flow_cont_qc')
            x_sites = fstats[fstats['Tot data yrs'] >= min_yrs].index.tolist()
        else:
            raise ValueError('Load some flow data in!')
    elif x_mtype == 'river_flow_disc_qc':
        ### gauging flow removal
        if 'river_flow_disc_qc' in self.mtypes:
            x_count = self.sel_ts('river_flow_disc_qc', x).groupby(level='site').count()
            x_sites = x_count[x_count >= min_obs].index.tolist()
        else:
            raise ValueError('Load some flow data in!')

    if y_mtype == 'river_flow_cont_qc':
        ### recorder flow removal
        if 'river_flow_cont_qc' in self.mtypes:
            sites1 = list(self.mtypes_sites['river_flow_cont_qc'])
            y_sites = Series(y)[Series(y).isin(sites1)].tolist()
        else:
            raise ValueError('Load some flow data in!')
    elif y_mtype == 'river_flow_disc_qc':
        ### gauging flow removal
        if 'river_flow_disc_qc' in self.mtypes:
            y_count = self.sel_ts('river_flow_disc_qc', y).groupby(level='site').count()
            y_sites = y_count[y_count >= min_obs].index.tolist()
        else:
            raise ValueError('Load some flow data in!')

    #### Create new hydro object with updated sites
    y_set = self.sel(mtypes=y_mtype, sites=y_sites)
    x_set = self.sel(mtypes=x_mtype, sites=x_sites)
    hydro1 = x_set.combine(y_set)

    #### Extract/create comparison dictionary
    if x is None:
        if buffer_dis is None:
            raise ValueError('If x is None, then buffer_dis should be defined.')
        hydro2 = hydro1._comp_by_buffer(buffer_dis)
        comp_dict = {i: list(hydro2.comp_dict[i]) for i in y_sites}
    elif isinstance(x, list):
        hydro2 = hydro1.copy()
        comp_dict = {i: x for i in y_sites}

    #### Run regression
    new_ts = DataFrame()
    reg = DataFrame()
    for j in comp_dict:
        x_data = hydro2.sel_ts(mtypes=x_mtype, sites=comp_dict[j], pivot=True)
        y_data = hydro2.sel_ts(mtypes=y_mtype, sites=j, pivot=True)
        y_data.name = j
        y_data = DataFrame(y_data)
        try:
            reg1, y_new1 = flow_reg_fun(x_data, y_data, min_obs=min_obs, p_val=p_val, logs=logs, make_ts=True, below_median=below_median)
            y_new1.name = j
            reg = concat([reg, reg1])
            new_ts = concat([new_ts, y_new1], axis=1)
        except:
            print('Site ' + str(j) + ' did not work with the regression.')
            pass
    new1 = self.add_data(new_ts, mtypes='river_flow_cont_qc', dformat='wide', add=False)
    return(new1, reg)



#######################################
#### Documentation for the functions

malf7d.__doc__ = malf_fun.__doc__






