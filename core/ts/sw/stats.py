# -*- coding: utf-8 -*-
"""
Flow stats functions.
"""
from os import path
from core.misc.misc import df_first_valid, df_last_valid
# from core.ecan_io.flow import rd_ts, rd_hydstra_csv
from core.ts.ts import tsreg
from scipy.signal import argrelmax
from core.stats.reg import lin_reg
from pandas import ExcelWriter, concat, Series, DataFrame, to_datetime, DateOffset
from core.allo_use import flow_ros, ros_freq
from numpy import in1d, nan, log, exp, argmax, max, where, mean, nanpercentile, round, floor


def flow_stats(x, below_median=False, export_path=None):
    """
    Function to run summary stats on time series flow data.

    Parameters
    ----------
    x : DataFrame or Series
        Pandas DataFrame or Series with a daily DateTimeIndex.
    below_median : bool
        Should the average number of days below the median be added?
    export_path : str
        Path where the results will be exported.

    Returns
    -------
    DataFrame
    """

    # Make sure the object is a data frame
    df_temp = DataFrame(x)
    df = df_temp.dropna(axis=1, how='all')

    # Run the stats
    zero_count = df[df == 0].count()
    d_count = df.count()
    min1 = df.min()
    max1 = df.max()
    mean1 = df.mean()
    max1 = df.max()
    median1 = df.median()
    quarter1 = df.apply(lambda x: nanpercentile(x, 25))
    quarter3 = df.apply(lambda x: nanpercentile(x, 75))
    max_time = df.idxmax()
    min_time = df.idxmin()
    first_time = df_first_valid(df)
    last_time = df_last_valid(df)
    periods = (last_time - first_time).astype('timedelta64[D]').astype('int') + 1
    mis_days = periods - d_count
    mis_days_ratio = mis_days / periods
    #    df2 = df[first_time.astype('str').values[0]:last_time.astype('str').values[0]]
    years = d_count / 365.25

    # prepare output
    #    df_names = df.columns.values
    row_names = ["Min", '25%', "Median", "Mean", '75%', "Max", "Start Date",
                 "End date", "Min date", "Max Date", "Zero days", "Missing days",
                 "Missing ratio", "Tot data yrs"]
    temp1 = concat([min1, quarter1, median1, mean1, quarter3, max1], axis=1).round(3)
    temp2 = concat(
        [first_time.astype('string'), last_time.astype('string'), min_time.astype('string'), max_time.astype('string'),
         zero_count, mis_days, mis_days_ratio.round(3), years.round(1)], axis=1)
    df2 = concat([temp1, temp2], axis=1)

    ### Add in the average number of days below median if specified
    if below_median:
        median_avg = median1.copy()
        row_names.extend(["Avg_days_below_median"])
        for i in median_avg.index:
            val = float(median_avg[i])
            t1 = df[i][df[i] <= val].count() / years[i]
            median_avg.loc[i] = t1.round(1)
        df2 = concat([df2, median_avg], axis=1)
    df2.columns = row_names
    df2.index.name = 'site'

    # Export data and return dataframe
    if isinstance(export_path, str):
        df2.to_csv(export_path)
    return (df2)


def malf7d(x, w_month='JUN', max_missing=90, malf_min=0.9, intervals=[10, 20, 30, 40], return_alfs=False,
           num_years=False, export_path=None, export_name_malf='malf.csv', export_name_alf='alf.csv',
           export_name_mis='alf_missing_data.csv'):
    """
    Function to create a 7 day mean annual low flow estimate from flow time
    series.

    Parameters
    ----------
    x : DataFrame or Series
        Pandas DataFrame or Series with a daily DateTimeIndex.
    w_month : str
        The month to start the water year in three upper case letters.
    max_missing : int
        The allowed missing data in a year for the alf7d calc.
    malf_min : float
        The minimum ratio of ALF data years to total years to calculate the MALF.
    intervals : list of int
        The intervals to calculate MALF over.
    return_alf : bool
        Should the ALFs and the number of missing days per year be returned in addition to the MALF?
    export_path : str
        The base path for the export files that will be saved.

    Returns
    -------
    DataFrame(s)
        When return_alfs is False, then the output is only a dataframe of the MALFs by intervals. When return_alfs is True, then the output will include the MALFs, the annual ALFs, the number of missing days per year, the number of days (out of 20) that have data surrounding the minimum flow value for that year, and the dates of the minimum flow per year.
    """

    mon_day_dict = {'JUN': 182, 'JAN': 1, 'FEB': 32, 'MAR': 61, 'APR': 92, 'MAY': 122, 'JUL': 153, 'AUG': 214,
                    'SEP': 245, 'OCT': 275, 'NOV': 306, 'DEC': 336}

    def malf_fun(df, intervals):

        malfs = []
        last_yr = df[-1:].index[0]
        for i in intervals:
            first_yr = last_yr - DateOffset(years=i)
            df10 = df[first_yr:]
            mis10 = floor(i * malf_min) <= df10.count()
            if mis10:
                malfs.extend([round(mean(df10), 3)])
            else:
                malfs.extend([nan])
        malfs.extend([round(mean(df), 3)])

        return (malfs)

    #    def day_june(df, dayofyear=182):
    #        day1 = df.dt.dayofyear
    #        over1 = day1 >= dayofyear
    #        under1 = day1 < dayofyear
    #        day2 = day1.copy()
    #        day2.loc[over1] = day1.loc[over1] - dayofyear
    #        day2.loc[under1] = 365 - dayofyear + day1.loc[under1]
    #        return(day2)

    ### Make sure the object is a data frame and regular
    df_temp = tsreg(DataFrame(x))
    df = df_temp.dropna(axis=1, how='all')
    df.columns = df.columns.astype(int)

    ### Rolling mean
    df2 = df.rolling(7, center=True).mean()

    ## count the number of days with data
    df2_nans = df2.rolling(20, center=True).count()

    ### Calc and filter ALFs
    n_days_yr = df2.fillna(0).resample('A-' + w_month).count()
    day_filter = n_days_yr.iloc[:, 0] >= 275
    n_days_yr2 = n_days_yr[day_filter]
    df3_res = df2.resample('A-' + w_month)
    df4 = df3_res.min()[day_filter]
    df_count1 = df3_res.count()[day_filter]
    df_missing = n_days_yr2 - df_count1
    df4[df_missing > max_missing] = nan
    df_count2 = df4.count()

    ### Find the minimum in each water year
    min_date = df2.resample('A-' + w_month).apply(lambda x: to_datetime(x.idxmin()))
    min_day = min_date.apply(lambda x: x.dt.dayofyear)
    mon_day = mon_day_dict[w_month]

    ## Determine if there are missing values near the min flow value
    n_days_min = min_date.copy()
    for i in df2_nans:
        index1 = min_date.loc[min_date[i].notnull(), i]
        val1 = df2_nans.loc[index1, i].resample('A-' + w_month).min()
        val1.name = i
        n_days_min.loc[:, i] = val1

    mis_min_bool = any(n_days_min < 13)

    if mis_min_bool:
        print(
        'Warning 1 - Some sites have significant amounts of missing values near the ALF! Check the DataFrame output for further info.')

    ## determine the number of min days per year within a month of the water year break point
    count_min_day = min_day.apply(lambda x: (x > (mon_day - 30)) & (x < (mon_day + 30))).sum()
    ratio_min_day = (count_min_day / df_count2).round(2)
    ratio_min_day_bool = ratio_min_day >= 0.25

    if any(ratio_min_day_bool):
        sites_prob = ratio_min_day[ratio_min_day_bool].index.tolist()
        print('Warning 2 - Site(s) ' + str(
            sites_prob) + ' have a significant amount of ALFs that fall at the end/beginning of the water year.')

    #    mean_date = min_date.apply(day_june, dayofyear=mon_day_dict[w_month]).mean()

    ## MALF calc
    malf = df4.apply(lambda x: Series(malf_fun(x, intervals)), axis=0).transpose()
    malf_col_names = ["MALF7D " + str(i) + " yrs" for i in intervals]
    malf_col_names.extend(["MALF7D all yrs"])
    malf.columns = malf_col_names
    if num_years:
        malf['Avg_num_days'] = nan
        for i in malf.index:
            val = float(malf.ix[i, 4])
            t1 = df[i][df[i] <= val]
            t2 = round(t1.resample('A-' + w_month).count().mean(), 1)
            malf.loc[i, 'Avg_num_days'] = t2

    ## Export data and return dataframes
    if return_alfs:
        if isinstance(export_path, str):
            malf.to_csv(path.join(export_path, export_name_malf))
            df4.round(3).to_csv(path.join(export_path, export_name_alf))
            df_missing.to_csv(path.join(export_path, export_name_mis))

        return ([malf, df4.round(3), df_missing, n_days_min, min_date])
    else:
        if isinstance(export_path, str):
            malf.to_csv(path.join(export_path, export_name_malf))

        return (malf)


def fre_accrual(csv_or_df, min_filter=False, min_yrs=4, max_missing=120, fre_mult=3, w_month='JUN', day_gap=5):
   """
   Function to calculate the Fre3 of a flow time series.
   """

   # if type(csv_or_df) is str:
   #     df = rd_hydstra_csv(csv_or_df, min_filter, min_yrs)
   #     df = df.dropna(axis=1, how='all')
   # else:

   df = csv_or_df
   df = df.dropna(axis=1, how='all')
   sites = df.columns.values
   howmany = len(sites)

   output = DataFrame(nan, index=sites, columns=['fre_val', 'fre_calc', 'mean_accrual', 'n_yrs'])

   medians = df.median()
   fre_val = medians * fre_mult

   for i in range(howmany):
       site = sites[i]
       threshold = fre_val[site]
       dfsite = df[[site]]

       start_date = df_first_valid(dfsite)
       start = start_date.iloc[0]
       #        start = to_datetime('2006-08-31')

       end_date = df_last_valid(dfsite)
       end = end_date.iloc[0]

       find_events = dfsite[start: end]
       # find_events = dfsite

       findgaps = find_events.dropna()
       findgaps['start_gap'] = findgaps.index
       findgaps['end_gap'] = findgaps['start_gap'].shift(periods=-1) + DateOffset(-1)
       findgaps['diff'] = abs(findgaps['start_gap'] - findgaps['end_gap'])
       gaps = findgaps.loc[findgaps['diff'] > '30 days']

       partitions = len(gaps) + 1

       find_events['Exceed'] = where(find_events[site] >= threshold, 1, 0)  # find days where threshold is exceeded
       find_events['DaysOver'] = find_events['Exceed'].rolling(
           window=6).sum()  # on how many days of the last 6 has the threshold been exceeded

       find_events['DaysOver'] = where((find_events['Exceed'] == 1), find_events['DaysOver'].fillna(1),
                                       find_events['DaysOver'].fillna(
                                           0))  # fill NaN values in DaysOver with 1 if flow over threshold, 0 if not

       find_events['Event'] = where((find_events['Exceed'] == 1) & (find_events['DaysOver'] == 1), 1,
                                    0)  # mark as a new event if the threshold is exceeded, and the 5 previous days were below the threshold

       # find the number of days in the year that are in the period of record, make sure it is at least 335
       n_days_yr = find_events[site].fillna(0).resample('A-' + w_month).count()
       day_filter = n_days_yr >= 335
       n_days_yr2 = n_days_yr[day_filter]

       events_year = find_events['Event'].resample('A-' + w_month).sum()[
           day_filter]  # find number of events that occurred each year

       data_days = find_events[site].resample('A-' + w_month).count()[day_filter]  # find number of days with data
       missing_days = n_days_yr2 - data_days  # find number of missing days per year
       events_year[
           missing_days > max_missing] = nan  # remove years with more than the maximum alllowance number of missing days
       n_years = events_year.count()  # find the length of the record

       # Fre calc
       fre_calc = events_year.mean()

       # Accrual periods

       if partitions == 1:
           if sum(find_events.Exceed) > 0:
               # calculate accrual period for the whole record
               find_accrual = find_events.loc[find_events['Exceed'] == 1]
               find_accrual['start_accrual'] = find_accrual.index
               find_accrual['end_accrual'] = find_accrual['start_accrual'].shift(periods=-1) + DateOffset(-1)
               find_accrual['accrual_time'] = abs(find_accrual['start_accrual'] - find_accrual['end_accrual'])
               accrual_periods = find_accrual.loc[find_accrual['accrual_time'] > (str(day_gap) + ' days')]
               mean_accrual = accrual_periods['accrual_time'].astype('timedelta64[D]').astype('int').mean()
       else:
           # split data into partitions if there are big gaps
           starts = gaps['end_gap'].tolist()
           starts.insert(0, start)
           ends = gaps['start_gap'].tolist()
           ends.insert(len(ends) + 1, end)
           storeaccrual = []

           for f in range(partitions):
               find_events_part = find_events[starts[f]:ends[f]]
               if sum(find_events_part.Exceed) > 0:
                   find_accrual = find_events_part.loc[find_events_part['Exceed'] == 1]
                   find_accrual['start_accrual'] = find_accrual.index
                   find_accrual['end_accrual'] = find_accrual['start_accrual'].shift(periods=-1) + DateOffset(-1)
                   find_accrual['accrual_time'] = abs(find_accrual['start_accrual'] - find_accrual['end_accrual'])
                   accrual_periods = find_accrual.loc[find_accrual['accrual_time'] > '5 days']
                   accrual_time = (accrual_periods['accrual_time'].astype('timedelta64[D]').astype('int')).tolist()
                   storeaccrual = storeaccrual + accrual_time
                   mean_accrual = mean(storeaccrual)

       # Put into output df
       output.loc[site, :] = [round(fre_val[site], 3), fre_calc, round(mean_accrual, 3), n_years]

   return (output)


def gauge_proc(gauge_flow, day_win=3, min_gauge=15, export=False, export_path='gaugings_flow.csv'):
    """
    Function to process gaugings data to determine the gauging periods with a minimum number of gaugings. Returns the gauging flows for those periods.

    Arguments:\n
    gauge_flow -- The output of the rd_henry function with sites_by_col=True.\n
    day_win -- The window period to look for concurrent gaugings over.\n
    min_gauge -- The minimum number of gaugings required.
    """

    ### Make the time series regular
    flow = tsreg(gauge_flow, freq='D')

    ### rolling window count and max
    flow_win = flow.rolling(day_win, min_periods=1, center=True).mean()

    flow_win_count = flow_win.count(axis=1)
    #    count10 = flow_win_count[flow_win_count > min_gauge]

    ### Determine local max from a smoothed spline
    flow_spline = flow_win_count.resample('2H').interpolate(method='spline', order=3, s=0.2)
    loc_max = argrelmax(flow_spline.values, order=48)
    flow_spline2 = flow_spline.ix[loc_max]
    flow_spline3 = flow_spline2[flow_spline2 > min_gauge]

    index_max = flow_spline3.index.round('D')

    ### Get the count and the final data
    flow2 = flow_win[flow_win.index.isin(index_max)].dropna(how='all', axis=1)

    ### Remove sites with only zero flow
    with_flow = flow2.mean() > 0.01
    flow3 = flow2.ix[:, with_flow.values]
    count_tot = flow3.count(axis=1)

    if export:
        flow3.to_csv(export_path)

    return ([flow3, count_tot])


def flow_reg(x, y, min_obs=10, p_val=0.05, logs=False, make_ts=False, cut_y=False, below_median=False, export=False, export_ts='est_flow.csv', export_reg='est_reg.csv'):
    """
    Function to make a regression between recorders and gaugings and optionally create a time series for the gaugings.

    Arguments:\n
    min_obs -- Minimum number of regressions values required.\n
    p_val -- Minimum p-value for the regressions.\n
    logs -- Should the variables be logged?\n
    make_ts -- Should a new time series be made from the regressions?\n
    cut -- Should the values above the max value from the regression be removed?
    """

    #### Filter out data that shouldn't be correlated
    x1 = x.replace(0, nan)
    y1 = y.replace(0, nan)

    #### Remove values below the median if desired
    if below_median:
        median1 = x1.median()

        new_df = DataFrame()
        for i in x1:
            med1 = median1.loc[i]
            grp = x1.loc[:, i]
            t1 = grp.loc[grp <= med1]
            new_df = concat([new_df, t1], axis=1)
        x2 = new_df
    else:
        x2 = x1.copy()

    #### Do the regressions
    t1 = lin_reg(x2, y1, log_x=logs, log_y=logs)
    t2 = [i[(i['n-obs'] >= min_obs)] for i in t1]
    t3 = [i for i in t2 if not i.empty]

    site_reg = []
    for g in y1.columns:
        site_reg1 = DataFrame((i[i['Y_loc'].values == g].values[0] for i in t3 if sum(i['Y_loc'].values == g)),
                              columns=t3[0].columns)
        site_reg.append(site_reg1)

    top1 = []
    for f in site_reg:
        p_index = f['p-value'] < p_val
        if sum(p_index) > 0:
            qual = argmax(((f.R2 - f.MANE))[p_index])
            best1 = f.iloc[qual, :]
            top1.append(best1)
    top = concat(top1, axis=1).T
    top.index = top['Y_loc']

    ### Create new time series from regressions
    if make_ts:
        y_new1 = []
        for i in top.index:
            reg1 = top.loc[i, :]
            if isinstance(x1, DataFrame):
                flow_x = x1[reg1['X_loc']]
            else:
                flow_x = x1
            flow_y = y1[reg1['Y_loc']]
            xy = concat([flow_x, flow_y], axis=1).dropna()

            ## Create new series
            if logs:
                y_new = exp(float(reg1.Slope) * log(flow_x) + float(reg1.Intercept))
            else:
                y_new = float(reg1.Slope) * flow_x + float(reg1.Intercept)

            ## remove negatives
            y_new[y_new < 0] = 0

            ## Remove all values above 'max_y'
            if cut_y:
                max_y = max(xy[reg1['Y_loc']]).astype('float32')
                y_new[y_new > (max_y + max_y * 0.2)] = nan

            ## Put back the original data
            y_new.loc[xy[reg1['Y_loc']].index] = y[i]

            ## Put into list container
            y_new.name = reg1['Y_loc']
            y_new1.append(y_new.round(3))

        ## Make into dataframe
        y_new2 = concat(y_new1, axis=1)

        ## Export if desired
        if export:
            y_new2.to_csv(export_ts)
            top.to_csv(export_reg)
        return ([top, y_new2])

    #### Export and return
    if export:
        top.to_csv(export_reg)
    return (top)


def est_low_flows(y, x, comp_dict, cut=True):
    """
    Estimate flows from regressions from two sites with a dictionary of the x: [y, log?].
    """

    x.columns = x.columns.astype(int)
    y.columns = y.columns.astype(int)
    df1 = DataFrame(nan, index=x.index, columns=[i for i in comp_dict])

    for i in comp_dict:
        y_site = int(i)
        x_site, logs = comp_dict[i]
        x_site = int(x_site)
        xy = concat([x[x_site], y[y_site]], axis=1)
        xy[xy <= 0] = nan
        xy = xy.dropna()
        max_y = max(xy[y_site])

        ## run regression and create new series
        reg1 = lin_reg(xy[x_site], xy[y_site], log_x=logs, log_y=logs)[0]
        if logs:
            y_new = exp(float(reg1.Slope) * log(x[x_site]) + float(reg1.Intercept))
        else:
            y_new = float(reg1.Slope) * x[x_site] + float(reg1.Intercept)

        ## remove negatives
        y_new[y_new < 0] = 0

        ## Remove all values about 'max_y'
        if cut:
            y_new[y_new > max_y] = nan

        ## Put back original data
        index1 = y[y_site].dropna().index
        y_new[index1] = y[y_site].dropna()

        df1[i] = y_new.round(3)

    return (df1)


def est_low_flows_reg(y, x, comp_dict, cut=True):
    """
    Estimate flows from regressions from two sites with a dictionary of the x: [y, log?].
    """

    x.columns = x.columns.astype(int)
    y.columns = y.columns.astype(int)
    df1 = []

    for i in comp_dict:
        y_site = int(i)
        x_site, logs = comp_dict[i]
        x_site = int(x_site)
        xy = concat([x[x_site], y[y_site]], axis=1)
        xy[xy <= 0] = nan
        xy = xy.dropna()

        ## run regression and create new series
        reg1 = lin_reg(xy[x_site], xy[y_site], log_x=logs, log_y=logs)[0]

        df1.append(reg1)

    ## Concat
    df2 = concat(df1)

    return (df2)


def flow_dur(flow, plot=False):
    """
    Function to estimate the flow duration of flow time series and optionally plot them.
    """

    flow1 = flow.dropna(how='all')
    rank1 = flow1.rank(axis=0, pct=True, ascending=False)

    fdc_sorted = []

    for i in flow1.columns:
        both2 = concat([flow1[i], rank1[i]], axis=1).dropna()
        both2.columns = ['flow1', 'rank1']
        both3 = both2.sort_values('rank1')
        fdc_sorted.extend([both3])

        if plot:
            both2.plot(x='rank1', y='flow1', kind='scatter', xlim=[0, 1], ylim=[0, max(both2.flow1)])

    return (fdc_sorted)


def summ_stats(df_lst, df_names, excel_file, below_median=True, num_years=True,
               flow_csv='S:/Surface Water/shared/base_data/flow/all_flow_data.csv'):
    """
    Summary flow stats function. Exports to an excel file. Input must be a list of time series dataframe(s). df_names must be a list of names associated with each dataframe (will be column headers).
    """

    stats_export = 'stats'
    ros_means_export = 'ros_means'
    ros_partial_export = 'ros_partial'
    ros_full_export = 'ros_full'
    ros_partial_norm_export = 'ros_partial_norm'
    ros_full_norm_export = 'ros_full_norm'

    ## Load all flow data
    all_flow = rd_ts(flow_csv)
    all_flow.columns = all_flow.columns.astype(int)

    #### Run loop for all three sets
    excel = ExcelWriter(excel_file)
    stats_out = []
    ros_means_out = []

    for i in range(len(df_lst)):
        flow = df_lst[i]

        ## Remove main sites from all flow data
        all_flow2 = all_flow.loc[:, ~in1d(all_flow.columns, flow.columns)]
        all_flow3 = concat([flow, all_flow2], axis=1, join='inner')

        ## Run stats for Huts and SH1 on the Pareora

        stats1 = flow_stats(flow, below_median=below_median)
        malf, alf, alf_mising = malf7d(flow, num_years=num_years)
        fre3 = fre_accrual(flow)

        # fdc, fdc = flow_dur(flow)

        stats2 = concat([stats1, malf.drop('MALF7D 40 yrs', axis=1), fre3], axis=1)
        stats_out.extend([stats2])

        #######################################
        #### Reliability of supply

        ros = flow_ros(flow.columns.values, flow_csv=all_flow3)
        ros_partial, ros_full = ros_freq(ros, period='summer')
        ros_partial_norm, ros_full_norm = ros_freq(ros, period='summer', norm=True)

        ros_partial_mean = ros_partial.mean().round(1)
        ros_full_mean = ros_full.mean().round(1)
        ros_partial_norm_mean = ros_partial_norm.mean().round(3)
        ros_full_norm_mean = ros_full_norm.mean().round(3)

        ros_means = concat([ros_partial_mean, ros_partial_norm_mean, ros_full_mean, ros_full_norm_mean], axis=1)
        ros_means.columns = ['ros_partial', 'ros_partial_norm', 'ros_full', 'ros_full_norm']

        ros_means_out.extend([ros_means])

        ######################################
        #### Export results
        out1 = df_names[i]

        stats2.to_excel(excel, out1 + '_' + stats_export)
        ros_means.to_excel(excel, out1 + '_' + ros_means_export)

        ros_partial.to_excel(excel, out1 + '_' + ros_partial_export)
        ros_full.to_excel(excel, out1 + '_' + ros_full_export)
        ros_partial_norm.to_excel(excel, out1 + '_' + ros_partial_norm_export)
        ros_full_norm.to_excel(excel, out1 + '_' + ros_full_norm_export)

    excel.save()
    return ([stats_out, ros_means_out])


def malf_reg(x, y, min_yrs=10, min_obs=10, w_month='JUN', max_missing=120, malf_min=0.9, intervals=[10, 20, 30, 40]):
    """
    A function to specifically correlate flows to estimate MALFs.
    """

    ### Only use x flow series with min_yrs
    m1 = malf7d(x, w_month=w_month, max_missing=max_missing, malf_min=malf_min, intervals=[min_yrs])
    x1 = x.loc[:, m1.iloc[:, 0].notnull()]

    ### Rolling 7D mean on x and y
    x_roll = x1.rolling(7, center=True).mean()
    y_roll = y.rolling(7, center=True).mean()

    ### Run the regressions
    reg1, df2 = flow_reg(x_roll, y_roll, min_obs=min_obs, make_ts=True, logs=False)

    ### Calc the 7DMALF
    def malf_fun(df, intervals):
        from numpy import mean, round, floor, nan

        malfs = []
        for i in intervals:
            df10 = df[-i:]
            mis10 = floor(len(df10) * round(1.0 - malf_min, 2)) >= sum(df10.isnull())
            if mis10:
                malfs.extend([round(mean(df10), 3)])
            else:
                malfs.extend([nan])
        malfs.extend([round(mean(df), 3)])

        return (malfs)

    ## Find the minimum in each water year and count the days
    n_days_yr = df2.fillna(0).resample('A-' + w_month).count()
    day_filter = n_days_yr.iloc[:, 0] >= 335
    n_days_yr2 = n_days_yr[day_filter]
    df3_res = df2.resample('A-' + w_month)
    df4 = df3_res.min()[day_filter]
    df_count1 = df3_res.count()[day_filter]
    df_missing = n_days_yr2 - df_count1
    df4[df_missing > max_missing] = nan
    #    df_count2 = df4.count()

    ## MALF calc
    malf = df4.apply(lambda x: Series(malf_fun(x, intervals)), axis=0).transpose()
    malf_col_names = ["MALF7D " + str(i) + " yrs" for i in intervals]
    malf_col_names.extend(["MALF7D all yrs"])
    malf.columns = malf_col_names

    ### Return results
    return ([reg1, malf])
