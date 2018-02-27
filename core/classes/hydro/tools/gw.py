# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:32:41 2017

@author: MichaelEK
"""

from pandas import DataFrame, concat, Series
from core.ts.sw.stats import flow_reg as flow_reg_fun
from numpy import in1d


def gwl_reg(self, y, x=None, y_mtype='aq_wl_disc_qc', x_mtype='aq_wl_cont_qc', buffer_dis=None, min_yrs=10, min_obs=10, logs=False, p_val=0.05, below_median=False):
    """
    Function to do a simple linear regression between manual GW sites and gwl recorder sites.
    The output is a hydro class object and the regression info.
    """
    #### Initial check
    if not isinstance(y, list):
        raise ValueError('y must be a list of site numbers')

    #### Remove data based on restriction parameters

    if x_mtype == 'aq_wl_cont_qc':
        ### recorder gwl removal
        if 'aq_wl_cont_qc' in self.mtypes:
            fstats = self.stats(mtypes='aq_wl_cont_qc')
            x_sites = fstats[fstats['Tot data yrs'] >= min_yrs].index.tolist()
        else:
            raise ValueError('Load some gwl data in!')
    elif x_mtype == 'aq_wl_disc_qc':
        ### gauging gwl removal
        if 'aq_wl_disc_qc' in self.mtypes:
            x_count = self.sel_ts('aq_wl_disc_qc', x).groupby(level='site').count()
            x_sites = x_count[x_count >= min_obs].index.tolist()
        else:
            raise ValueError('Load some gwl data in!')

    if y_mtype == 'aq_wl_cont_qc':
        ### recorder gwl removal
        if 'aq_wl_cont_qc' in self.mtypes:
            sites1 = list(self.mtypes_sites['aq_wl_cont_qc'])
            y_sites = Series(y)[Series(y).isin(sites1)].tolist()
        else:
            raise ValueError('Load some gwl data in!')
    elif y_mtype == 'aq_wl_disc_qc':
        ### gauging gwl removal
        if 'aq_wl_disc_qc' in self.mtypes:
            y_count = self.sel_ts('aq_wl_disc_qc', y).groupby(level='site').count()
            y_sites = y_count[y_count >= min_obs].index.tolist()
        else:
            raise ValueError('Load some gwl data in!')

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
    new1 = self.add_data(new_ts, mtypes='aq_wl_cont_qc', dformat='wide', add=False)
    return(new1, reg)
