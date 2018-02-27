# -*- coding: utf-8 -*-
"""
General time series tools.
"""
#from pandas.tseries.resample import Resampler
#from pandas.core.groupby import GroupBy


def resample(self, resample_code='A-JUN', fun='mean'):
    """
    Time series resampling function. Returns a Hydro class object with resampled data.

    Parameters
    ----------
    resample_code : str
        The Pandas resampling code for the resampling process (e.g. 'A' for annual, 'A-JUN' for annual ending in June (water year), 'D' for day, 'W' for week, 'M' for month')
    fun : str
        The function that should be applied. Any function that Pandas can handle in a groupby object.
    """
    from pandas.core.groupby import SeriesGroupBy, GroupBy
    from pandas import Grouper

    ### Set up parameters
    if fun in GroupBy.__dict__.keys():
        fun1 = GroupBy.__dict__[fun]
    elif fun in SeriesGroupBy.__dict__.keys():
        fun1 = SeriesGroupBy.__dict__[fun]
    else:
        raise ValueError('Need to use the right function name.')

    ### Run resampling
    data = self.data

    df1 = data.groupby([Grouper(level='mtype'), Grouper(level='site'), Grouper(level='time', freq=resample_code)])
    df2 = fun1(df1)

    ### Recreate hydro class object
    new1 = self.add_data(df2.reset_index(), 'time', 'site', 'mtype', 'data', dformat='long', add=False)
    return(new1)


def stats(self, mtypes, sites=None, below_median=False):
    """
    Function to produce stats for specific mytpes.

    mtypes : str
        A single str easurement type.
    sites : str, int, or list
        A str, int, or list of sites that should be included.
    below_median : bool
        Specific for the 'flow' stats.
    """
    from core.ts.sw.stats import flow_stats
    from core.ts.met.met_stats import precip_stats

    if isinstance(mtypes, str):
        if 'river_flow_cont' in mtypes:
            data = self.sel_ts(mtypes=['river_flow_cont'], sites=sites, pivot=True)
            if data.index.inferred_freq != 'D':
                data = data.resample('D').mean()
            stats1 = flow_stats(data, below_median=below_median)
            return(stats1)
        if 'atmos_precip_cont' in mtypes:
            data = self.sel_ts(mtypes=['atmos_precip_cont'], sites=sites)
            data.index = data.index.droplevel('mtype')
            stats1 = precip_stats(data)
            return(stats1)
        if 'aq_wl_cont' in mtypes:
            data = self.sel_ts(mtypes=['aq_wl_cont'], sites=sites)
            data.index = data.index.droplevel('mtype')
            stats1 = precip_stats(data)
            return(stats1)
    else:
        raise ValueError('No stats for the mtype')


