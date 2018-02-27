# -*- coding: utf-8 -*-
"""
General time series functions.
"""
from pandas.core.groupby import SeriesGroupBy, GroupBy
from pandas import infer_freq, TimeGrouper, Timestamp, Grouper, DataFrame
from numpy import nan, array
from os import path
from core.misc.misc import time_switch


def pd_grouby_fun(fun_name):
    """
    Function to make a function specifically to be used on pandas groupby objects from a string code of the associated function.
    """

    if fun_name in GroupBy.__dict__.keys():
        fun1 = GroupBy.__dict__[fun_name]
    elif fun_name in SeriesGroupBy.__dict__.keys():
        fun1 = SeriesGroupBy.__dict__[fun_name]
    else:
        raise ValueError('Need to use the right function name.')
    return (fun1)


def pd_ts_code(period, n_periods):
    """
    Function to convert common time period names and an associated number of periods to a pandas time string.

    period -- A string of common time period names (e.g. 'day', 'month', 'year', 'water year').\n
    n_periods -- The number of periods.
    """

    time_code = str(n_periods) + time_switch(period)

    return (time_code)


def tsreg(ts, freq=None, interp=False):
    """
    Function to regularize a time series object (pandas).
    The first three indeces must be regular for freq=None!!!

    ts -- pandas time series dataframe.\n
    freq -- Either specify the known frequency of the data or use None and
    determine the frequency from the first three indices.\n
    interp -- Should linear interpolation be applied on all missing data?
    """

    if freq is None:
        freq = infer_freq(ts.index[:3])
    ts1 = ts.resample(freq).mean()
    if interp:
        ts1 = ts1.interpolate('time')

    return (ts1)


def w_agg(x, fun='sum', axis=1):
    """
    Aggregates a dataframe of data into a single series.
    """
    if x.index.dtype_str != 'datetime64[ns]':
        print("Not time series!")
    if fun == 'mean':
        agg1 = x.mean(axis=axis)
    elif fun == 'sum':
        agg1 = x.sum(axis=axis)

    return (agg1)


def grp_ts_agg(df, grp_col, ts_col, freq_code):
    """
    Simple function to aggregate time series with dataframes with a single column of sites and a column of times.

    Arguments:\n
    df -- dataframe with a datetime column.\n
    grp_col -- Column name that contains the sites.\n
    ts_col -- The column name of the datetime column.\n
    freq_code -- The pandas frequency code for the aggregation (e.g. 'M', 'A-JUN').\n
    agg_fun -- Either 'mean' or 'sum'.
    """

    df1 = df.copy()
    if type(df[ts_col].iloc[0]) is Timestamp:
        df1.set_index(ts_col, inplace=True)
        if type(grp_col) is list:
            grp_col.extend([TimeGrouper(freq_code)])
        else:
            grp_col = [grp_col, TimeGrouper(freq_code)]
        df_grp = df1.groupby(grp_col)
        return (df_grp)
    else:
        print('Make one column a timeseries!')


def interp_resample(df, res_code, level=None):
    """
    Function to properly set up a resampling class for discrete data. This assumes a linear interpolation between data points.

    Parameters
    ----------
    df: DataFrame
        DataFrame with a time index.
    res_code: str
        Pandas resampling code. e.g. 'D'.
    level: str or int
        For a MultiIndex, level (name or number) to use for resampling. Level must be datetime-like.

    Returns
    -------
    Resampling class
    """
    df1 = (df + df.shift(-1))/2
    out1 = df1.resample(res_code, level=level)
    return out1


def w_resample(x, period='water year', n_periods=1, fun='sum', min_ratio=0.75, agg=False, agg_fun='mean', digits=3, interp=False, discrete=False, export=False, export_path='', export_name='precip_stats.csv'):
    """
    NEED TO UPDATE! CLEAN AND INCLUDE THE discrete PARAMETER USING interp_resample.

    Function to resample time series precip or flow data.

    Arguments:\n
    period -- The length of time that the data should be resampled. Options are: 'days', 'weeks', 'months', 'years', and 'water years' (and the singular).\n
    n_periods -- the number of 'periods' to aggregate over. E.g. when 'period' is set to 'months' and n_periods is equal to 2, then it aggregates over a two month period.\n
    fun -- The function that should be used to resample ('mean' or 'sum').\n
    min_ratio -- The minimum proportion of data to total days for the final estimate.\n
    agg -- Should the columns be aggregated together?\n
    agg_fun -- The function for the aggregation ('mean' or 'sum').\n
    interp -- Should missing days be interpolated by the period mean?\n
    export -- Should the results be exported?\n
    export_path -- Directory where the results will be exported.\n
    export_name -- The name of the csv file to be exported.
    """

    # Make sure the object is a data frame and regular
    df = x.copy()

    # Warn about duplicate stations
    if sum(df.columns.duplicated()) > 0:
        print("Duplicate stations!!! Please check!!!")
        df = df.loc[:, df.columns.duplicated() == False]

    # Aggregate stations if desired
    if agg:
        df = w_agg(df, fun=agg_fun)

    # Resample
    df_grp = df.resample(str(n_periods) + time_switch(period))

    # Determine number of days in each period
    count1 = df_grp.count()
    tot_count1 = df_grp.size()
    diff_count = count1.sub(tot_count1, axis='index').abs()
    max1 = max(tot_count1)
    tot_count1.iloc[[0, -1]] = max1
    ratio1 = count1.div(tot_count1, axis='index')

    # Apply function to resample
    if fun == 'mean':
        res1 = df_grp.mean()
        if isinstance(res1, DataFrame):
            for i in res1.columns:
                res1.loc[ratio1[i] < min_ratio, i] = nan
        else:
            res1.loc[ratio1 < min_ratio] = nan
    elif fun == 'sum':
        res1 = df_grp.sum()
        mean1 = df_grp.mean()
        if isinstance(res1, DataFrame):
            for i in res1.columns:
                if interp:
                    res1.loc[ratio1[i] < 0.95, i] = res1.loc[ratio1[i] < 0.95, i] + mean1.loc[ratio1[i] < 0.95, i] * \
                                                                                    diff_count.loc[ratio1[i] < 0.95, i]
                res1.loc[ratio1[i] < min_ratio, i] = nan
        else:
            if interp:
                res1.loc[ratio1 < 0.95] = res1.loc[ratio1 < 0.95] + mean1.loc[ratio1 < 0.95] * diff_count.loc[
                    ratio1 < 0.95]
            res1.loc[ratio1 < min_ratio] = nan

    # Export data and return dataframe
    if export:
        res1.round(digits).to_csv(path.join(export_path, export_name))
    return (res1.round(digits))


def res(ts, dformat, wide_index=[0, 1], period='water year', n_periods=1, fun='sum', min_ratio=0.75, agg=False,
        agg_fun='mean', digits=3, interp=False, export=False, export_path='', export_name='precip_stats.csv'):
    """
    Function to resample time series precip or flow data.

    Arguments:\n
    period -- The length of time that the data should be resampled. Options are: 'days', 'weeks', 'months', 'years', and 'water years' (and the singular).\n
    n_periods -- the number of 'periods' to aggregate over. E.g. when 'period' is set to 'months' and n_periods is equal to 2, then it aggregates over a two month period.\n
    fun -- The function that should be used to resample ('mean' or 'sum').\n
    min_ratio -- The minimum proportion of data to total days for the final estimate.\n
    agg -- Should the columns be aggregated together?\n
    agg_fun -- The function for the aggregation ('mean' or 'sum').\n
    interp -- Should missing days be interpolated by the period mean?\n
    export -- Should the results be exported?\n
    export_path -- Directory where the results will be exported.\n
    export_name -- The name of the csv file to be exported.
    """

    # Make sure the object is a data frame and regular
    df = ts.copy()

    if dformat == 'wide':

        # Warn about duplicate stations
        if sum(df.columns.duplicated()) > 0:
            raise ValueError("Duplicate stations!!! Please check!!!")

        # Aggregate stations if desired
        if agg:
            df = w_agg(df, fun=agg_fun)

        # Resample
        df_grp = df.resample(str(n_periods) + time_switch(period))

    elif dformat == 'long':

        # Warn about duplicate stations
        if sum(df.index.levels[0].duplicated()) > 0:
            raise ValueError("Duplicate stations!!! Please check!!!")

        df_grp = df.groupby(
            [Grouper(level=wide_index[0]), Grouper(level=wide_index[1], freq=str(n_periods) + time_switch(period))])

    else:
        raise ValueError('dformat must be either long or wide')

    # Determine number of days in each period
    count1 = df_grp.count()
    tot_count1 = df_grp.size()
    diff_count = count1.sub(tot_count1, axis='index').abs()
    max1 = max(tot_count1)
    tot_count1.iloc[[0, -1]] = max1
    ratio1 = count1.div(tot_count1, axis='index')

    # Apply function to resample
    if dformat == 'wide':

        if fun == 'mean':
            res1 = df_grp.mean()
            if isinstance(res1, DataFrame):
                for i in res1.columns:
                    res1.loc[ratio1[i] < min_ratio, i] = nan
            else:
                res1.loc[ratio1 < min_ratio] = nan
        elif fun == 'sum':
            res1 = df_grp.sum()
            mean1 = df_grp.mean()
            if isinstance(res1, DataFrame):
                for i in res1.columns:
                    if interp:
                        res1.loc[ratio1[i] < 0.95, i] = res1.loc[ratio1[i] < 0.95, i] + mean1.loc[ratio1[i] < 0.95, i] * \
                                                                                        diff_count.loc[
                                                                                            ratio1[i] < 0.95, i]
                    res1.loc[ratio1[i] < min_ratio, i] = nan
            else:
                if interp:
                    res1.loc[ratio1 < 0.95] = res1.loc[ratio1 < 0.95] + mean1.loc[ratio1 < 0.95] * diff_count.loc[
                        ratio1 < 0.95]
                res1.loc[ratio1 < min_ratio] = nan

    elif dformat == 'long':

        if fun == 'mean':
            res1 = df_grp.mean()
        if fun == 'sum':
            res1 = df_grp.sum()
        res1[ratio1 < min_ratio] = nan

    # Export data and return dataframe
    if export:
        res1.round(digits).to_csv(path.join(export_path, export_name))
    return (res1.round(digits))


def ts_comp(df, freq='water year', n_periods=1, fun='mean', min_ratio=0.75, digits=3, start='1900', end='2020'):
    """
    Function to compare two sets of time series over a given period.
    """

    ### Resample data
    res1 = w_resample(df, freq, n_periods, fun, min_ratio, digits=digits)[start:end]

    ### Run stats
    res1_count = res1.count()
    res1_mean = res1.mean()
    res1_mean_inner = res1.dropna().mean()

    ### Organize results
    df2 = DataFrame(array([res1_count, res1_mean, res1_mean_inner]),
                    index=['n_yrs', 'Mean_long_record', 'Mean_short_record'], columns=res1.columns)

    return (df2)
