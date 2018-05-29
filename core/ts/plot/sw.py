# -*- coding: utf-8 -*-
"""
Surface water plotting functions.
"""


def hydrograph_plot(flow, precip=None, x_period='day', x_n_periods=1, time_format='%d-%m %H', x_rot=45, alpha=0.6, flow_units='m^{3}/s', precip_units='mm', x_grid=False, export_path=None):
    """
    Plotting function for hydrographs with/without a precipitation data set.

    Arguments:\n
    flow -- Pandas Dataframe time series with sites as columns.\n
    precip -- Pandas Dataframe time series with sites as columns.\n
    x_period -- The gross period for the axis.\n
    x_n_period -- The number of intervals per x_period.\n
    time_format -- The Date-time format for the x-axis in ISO format.\n
    x_rot -- The rotation of the x-axis labels.\n
    alpha -- The transparency of the precip plots.\n
    flow_units -- The units that will be shown on the flow y-axis.\n
    precip_units -- The units that will be shown on the precip y-axis.\n
    x_grid -- Should a vertical grid be drawn?
    """
    import matplotlib.pyplot as plt
    from numpy import arange, in1d, where
    from seaborn import set_style, despine, set_context, color_palette
    from core.ts import w_resample

    set_style("white")
    set_style("ticks")
    set_context('poster')
    pal1 = color_palette()

    #### Set up index and x scales
    x_index = flow.index
    x_index2 = w_resample(flow, x_period, x_n_periods).index
    x_axis = arange(0, len(flow.index))
    end = x_axis[-1]
    xticks_num = where(in1d(x_index, x_index2))[0]
    xticks = x_index[in1d(x_index, x_index2)].strftime(time_format)

    #### Make the plots
    u1 = flow.reset_index().iloc[:, 1:]
    ax1 = u1.plot(x_axis, rot=x_rot)
    if precip is not None:
        ax2 = ax1.twinx()
        u2 = precip.reset_index().iloc[:, 1:]
        x_axis2 = x_axis[in1d(flow.index, precip.index)]
        precip_names = u2.columns
        bar_width = ((x_axis2[1] - x_axis2[0]) * 0.75)/len(precip_names)
        for i in arange(0, len(precip_names)):
            ax2.bar(x_axis2 + bar_width*i, u2.iloc[:,i], bar_width, color=pal1[i], label=precip_names[i], alpha=alpha)
        ax2.set_ylabel('Precipitation $(' + precip_units + ')$')
        ax2.legend(loc='upper right')
        despine(right=False)
    else:
        despine()

    #### Adjust the axes and legends
    ax1.xaxis.set_ticks(xticks_num)
    plt.xticks(xticks_num, xticks)
    plt.xlim(0, end)
    ax1.set_ylabel('Flow $(' + flow_units + ')$')
    ax1.legend(loc='upper left')
    if x_grid:
        ax1.xaxis.grid(True, linewidth=0.5, linestyle='--')
        despine(right=False, top=False)

    plt.tight_layout()
    plot = ax1.get_figure()

    #### Save plot
    if isinstance(export_path, str):
        plot.savefig(export_path)

    return(plot)


def mon_boxplot(ts, col='all', dtype='flow', fun='mean', min_yrs=10, export_path='', export_name='flow_boxplot.png'):
    """

    """
    from os import path
    from calendar import month_abbr
    import seaborn as sns
    import matplotlib.pyplot as plt
    from core.ts import w_resample
    from core.ts.sw import flow_stats
    from pandas import Series

    ylabel_dict = {'flow': 'Mean Flow $(m^{3}/s)$', 'precip': 'Precipitation $(mm/month)$', 'PET': 'Potential Evapotranspiration $(mm/month)$', 'diff': 'Precipitation - PET $(mm/month)$'}
    sns.set_style("whitegrid")
    sns.set_context('poster')
    col_dict = {'flow': sns.color_palette("pastel")[0], 'precip': sns.color_palette("pastel")[0], 'PET': sns.color_palette("pastel")[1], 'diff': sns.light_palette("blue")[0]}
    color = col_dict[dtype]

    if type(ts) is Series:
        ts2 = ts
    else:
        stats1 = flow_stats(ts)
        index1 = stats1.loc['Tot data yrs', :] >= min_yrs
        sites = ts.columns[index1]
        ts2 = ts[sites]
    tsm = w_resample(ts2, period='month', fun=fun)
    tsm['month_num'] = tsm.index.month
    tsm['month'] = tsm['month_num'].apply(lambda x: month_abbr[x])

    if type(ts) is Series:
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.boxplot(x='month', y=0, data=tsm, color=color, showfliers=False)
        plt.ylabel(ylabel_dict[dtype])
        plt.xlabel('Month')
        plt.tight_layout()
        plot2 = ax.get_figure()
        plot2.savefig(path.join(export_path, export_name))
    elif col is 'all':
        for i in sites:
            fig, ax = plt.subplots(figsize=(15, 10))
            sns.boxplot(x='month', y=i, data=tsm, color=color, showfliers=False)
            plt.ylabel(ylabel_dict[dtype])
            plt.xlabel('Month')
            plt.tight_layout()
            plot2 = ax.get_figure()
            plot2.savefig(path.join(export_path, str(i) + '_' + export_name))
    elif type(col) is list:
        for i in col:
            fig, ax = plt.subplots(figsize=(15, 10))
            sns.boxplot(x='month', y=i, color=color, data=tsm)
            plt.ylabel(ylabel_dict[dtype])
            plt.xlabel('Month')
            plt.tight_layout()
            plot2 = ax.get_figure()
            plot2.savefig(path.join(export_path, str(i) + '_' + export_name))


def dual_mon_boxplot(ts1, ts2, col='all', dtype='flow', fun='mean', min_yrs=10, export_path='', export_name='flow_dual_boxplot.png'):
    """

    """
    from calendar import month_abbr
    from os import path
    from pandas import concat, melt, DatetimeIndex
    import seaborn as sns
    import matplotlib.pyplot as plt
    from core.ts.ts import w_resample
    from core.ts.sw.stats import flow_stats

    ### Set base parameters
    ylabel_dict = {'flow': 'Mean Flow $(m^{3}/s)$', 'precip': 'Precipitation $(mm/month)$', 'ET': 'Evapotranspiration $(mm/month)$'}
    sns.set_style("whitegrid")
    sns.set_context('poster')
    cat = ['Recorded', 'Naturalised']

    ### Process data
    ts1 = ts1.dropna(axis=1, how='all')
    stats1 = flow_stats(ts1)
    index1 = stats1.loc[:, 'Tot data yrs'] >= min_yrs
    sites = ts1.columns[index1.values]
    ts1_sel = ts1[sites]
    ts2_sel = ts2[sites]
    ts1_selm = w_resample(ts1_sel, period='month', fun=fun)
    ts2_selm = w_resample(ts2_sel, period='month', fun=fun)
    ts_all = concat([ts1_selm, ts2_selm], axis=1, join='inner')

    ### Make the figures
    if col is 'all':
        for i in sites:
            ts_site = ts_all[i]
            ts_site.columns = cat
            ts_site2 = melt(ts_site.reset_index(), id_vars='Time', value_vars=cat, var_name='type')
            ts_site2['month_num'] = DatetimeIndex(ts_site2.Time).month
            ts_site2['month'] = ts_site2['month_num'].apply(lambda x: month_abbr[x])
            ts_site2 = ts_site2.sort_values(['month_num', 'type'], ascending=[True, False], axis='index')
            fig, ax = plt.subplots(figsize=(15, 10))
            sns.boxplot(x='month', y='value', hue='type', data=ts_site2, showfliers=False)
            plt.ylabel(ylabel_dict[dtype])
            plt.xlabel('Month')
            plt.legend()
            plt.tight_layout()
            plot2 = ax.get_figure()
            plot2.savefig(path.join(export_path, str(i) + '_' + export_name))
    elif type(col) is list:
        for i in col:
            ts_site = ts_all[i]
            ts_site.columns = cat
            ts_site2 = melt(ts_site.reset_index(), id_vars='Time', value_vars=cat, var_name='type')
            ts_site2['month_num'] = DatetimeIndex(ts_site2.Time).month
            ts_site2['month'] = ts_site2['month_num'].apply(lambda x: month_abbr[x])
            ts_site2 = ts_site2.sort_values(['month_num', 'type'], ascending=[True, False], axis='index')
            fig, ax = plt.subplots(figsize=(15, 10))
            sns.boxplot(x='month', y='value', hue='type', data=ts_site2, showfliers=False)
            plt.ylabel(ylabel_dict[dtype])
            plt.xlabel('Month')
            plt.legend()
            plt.tight_layout()
            plot2 = ax.get_figure()
            plot2.savefig(path.join(export_path, str(i) + '_' + export_name))


def multi_yr_barplot(ts1, ts2, col='all', names_dict=None, dtype='flow', fun='mean', single=False, start='1970', end='2015', alf=False, export_path='', export_name='flow_multi_barplot.png'):
    """

    """
    from os import path
    from numpy import repeat
    from pandas import concat, melt, DatetimeIndex
    import seaborn as sns
    import matplotlib.pyplot as plt
    from core.ts import w_resample, tsreg
    from core.ts.sw import flow_stats, malf7d

    ### Set base parameters
    ylabel_dict = {'flow': 'Mean Flow $(m^{3}/s)$', 'precip': 'Precipitation $(mm/year)$', 'ET': 'Evapotranspiration $(mm/year)$', 'both': '$mm/year$'}
    sns.set_style("whitegrid")
    sns.set_context('poster')
    cat_dict = {'flow': ['Recorded', 'Naturalised'], 'both': ['Precipitation', 'Potential Evapotranspiration']}
    cat = cat_dict[dtype]

    ### assign names if needed
    if type(names_dict) is dict:
        col_names = [names_dict[x] for x in col]

    ### Process data
    if col is not 'all':
        ts1_sel = ts1[col]
        ts2_sel = ts2[col]
        ts1_sel.columns = ts2_sel.columns = col_names
    else:
        ts1_sel = ts1
        ts2_sel = ts2

    ts1_yr = w_resample(ts1_sel, period='water year', fun=fun)
    ts2_yr = w_resample(ts2_sel, period='water year', fun=fun)

    if alf:
        malf1, ts1_yr, alf_mis1 = malf7d(ts1_sel)
        malf1, ts2_yr, alf_mis1 = malf7d(ts2_sel)

    ts_all = tsreg(concat([ts1_yr, ts2_yr], axis=1, join='inner', keys=cat).dropna(how='all'))
    ts_all = ts_all[start:end]

    if single:
        ts_non = ts_all.unstack().reset_index().drop('level_1', axis=1)
        time_col = ts_all.index.name
        ts_non['year'] = DatetimeIndex(ts_non[time_col]).year
        ts_non.columns = ['site', time_col, 'value', 'year']
        ts_non.site = ts_non.site.astype('str')

        fig, ax = plt.subplots(figsize=(15, 10))
        sns.barplot(x='year', y='value', hue='site', data=ts_non,   palette='pastel')
        lgd = plt.legend(loc='best')

    else:
        ts_non = ts_all[cat[0]].unstack().reset_index()
        time_col = ts_all.index.name
        ts_non['year'] = DatetimeIndex(ts_non[time_col]).year
        ts_non.columns = ['site', time_col, 'value', 'year']
        ts_non.site = ts_non.site.astype('str')

        ts_nat = ts_all[cat[1]].unstack().reset_index()
        ts_nat['year'] = DatetimeIndex(ts_nat[time_col]).year
        ts_nat.columns = ['site', time_col, 'value', 'year']
        ts_nat.site = ts_nat.site.astype('str')

        order1 = ts1_sel.columns.astype('str').tolist()

        ### Make the figures
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.barplot(x='year', y='value', hue='site', data=ts_nat,   palette='pastel', hue_order=order1)
        sns.barplot(x='year', y='value', hue='site', data=ts_non, saturation=0.75, hue_order=order1)
        plt.plot([0], marker='None', linestyle='None', label='dummy')
        handles, lbs = ax.get_legend_handles_labels()
        handles2 = [handles[0]] + handles[len(col)+1:] + handles[:len(col)+1]
        labels2 = [cat[0]] + order1 + [cat[1]] + order1
        lgd = plt.legend(handles2, labels2, loc='center left', bbox_to_anchor=(1, 0.5))

    xticks = ax.get_xticks()
    if len(xticks) > 15:
        for label in ax.get_xticklabels()[::2]:
            label.set_visible(False)
        ax.xaxis_date()
        fig.autofmt_xdate(ha='center')
#        plt.tight_layout()

#    ax.xaxis_date()
#    fig.autofmt_xdate(ha='center')
    if alf:
        plt.ylabel('7DALF $(m^{3}/s)$')
    else:
        plt.ylabel(ylabel_dict[dtype])
    plt.xlabel('Water Year')
#    lgd = plt.legend([handles[i] for i in order2], [r'$D_{etc}$'] + cat + [r'$A_{etc}$'] + cat, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
#    plt.tight_layout()
    plot2 = ax.get_figure()
    plot2.savefig(path.join(export_path, export_name), bbox_extra_artists=(lgd,), bbox_inches='tight')


def reg_plot(x, y, freq='day', n_periods=1, fun='mean', min_ratio=0.75, digits=3, x_max=None, y_max=None, logs=False, export=False, export_path='flow_reg.png'):
    """
    Function to plot regressions with confidence intervals.
    """
    from os import path
    from numpy import nan, log
    from core.stats.reg import lin_reg
    import seaborn as sns
    from core.ts.ts import w_resample
    from pandas import concat
    import matplotlib.pyplot as plt

    ### Resample and combine
    x1 = w_resample(x, freq, n_periods, fun, min_ratio, digits)
    y1 = w_resample(y, freq, n_periods, fun, min_ratio, digits)

    xy1 = concat([x1, y1], axis=1, join='inner').dropna()
    col_names = [str(x.columns[0]), str(y.columns[0])]
    xy1.columns = col_names
    if logs:
        xy1 = xy1.apply(log)

    ### Make regression
    if x_max and y_max is not None:
        xy2 = xy1.copy()
        xy2.ix[xy2.iloc[:, 0] > x_max, 0] = nan
        xy2.ix[xy2.iloc[:, 1] > y_max, 1] = nan
        xy2 = xy2.dropna()
        reg1 = lin_reg(xy2.iloc[:, 0], xy2.iloc[:, 1])[0]
    else:
        xy2 = xy1.copy()
        reg1 = lin_reg(xy2.iloc[:, 0], xy2.iloc[:, 1])[0]

    ### plot regression
    sns.set_style("whitegrid")
    sns.set_context('poster')
    fig, ax = plt.subplots(figsize=(15, 10))
#    if log_x:
#        ax.set(xscale="log")
#    if log_y:
#        ax.set(yscale="log")
    sns.regplot(x=col_names[0], y=col_names[1], data=xy2.reset_index(), truncate=False, fit_reg=True, ax=ax)
    x_lim = ax.xaxis.get_view_interval()
    y_lim = ax.yaxis.get_view_interval()
    if x_max is not None:
        ax.set_xlim(0, x_max)
#    else:
#        ax.set_xlim(0, x_lim[1])
    if y_max is not None:
        ax.set_ylim(0, y_max)
#    else:
#        ax.set_ylim(0, y_lim[1])
    plt.xlabel(col_names[0] + ' $(m^{3}/s)$')
    plt.ylabel(col_names[1] + ' $(m^{3}/s)$')

    plot2 = ax.get_figure()
    if export:
        plot2.savefig(path.join(export_path), bbox_inches='tight')

    return(plot2, reg1)













