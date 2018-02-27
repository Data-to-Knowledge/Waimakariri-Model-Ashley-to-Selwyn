# -*- coding: utf-8 -*-
"""
Misc functions for the hydro class.
"""
from pandas import DataFrame, Series, DatetimeIndex, to_datetime, MultiIndex, concat
from numpy import array, ndarray, in1d, unique, append, nan, argmax

### Internal checking functions

def _check_mtypes_sites(self, mtype):
    d1 = {i: all(in1d(mtype, self.sites_mtypes[i])) for i in self.sites_mtypes}
    return(d1)


### Base stats for the default view of the class (once data has been loaded)

def _base_stats_fun(self):
    grp1 = self.data.groupby(level=['mtype', 'site'])
    start = grp1.apply(lambda x: x.first_valid_index()[2])
    start.name = 'start_time'
    end = grp1.apply(lambda x: x.last_valid_index()[2])
    end.name = 'end_time'
    stats1 = grp1.describe()[['min', '25%', '50%', '75%', 'mean', 'max', 'count']].round(2)
    out1 = concat([stats1, start, end], axis=1)
    setattr(self, '_base_stats', out1)

### mtypes checking for imprting tools

def _mtype_check(self, name):
    cond = name in self.mtypes
    return(cond)


