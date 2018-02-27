# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 13:57:22 2016

@author: Kate and Mike
"""


def allo_band_plt(csv_path, river):

    from pandas import read_csv
    import matplotlib.pyplot as plt
    matplotlib.style.use('ggplot')

    Alloc = read_csv(csv_path, index_col=(0), header=[0,1], skiprows=[2])
    Alloc = Alloc.fillna(0)
    Alloc = Alloc.rename(columns={'S&D no min flow': 'Stock and domestic, no minimum flow'})

    #Irrigation use and allocation band plot
    Alloc1 = Alloc.sum(axis=1, level=1)

    AllocPlot1 = Alloc1.plot.bar(stacked=True,colormap='Set3')
    AllocPlot1 = plt.xlabel('Year')
    AllocPlot1 = plt.ylabel('Allocation (L/s)')
    AllocPlot1 = plt.title(river + ' allocation')
    AllocPlot2 = plt.legend(loc='upper left', title="Water use & allocation band")

    #Irrigation source plot
    Alloc2 = Alloc.sum(axis=1, level=0)

    AllocPlot2 = Alloc2.plot.area(stacked=True)
    AllocPlot2 = plt.xlabel('Year')
    AllocPlot2 = plt.ylabel('Allocation (L/s)')
    AllocPlot2 = plt.title(river + ' allocation')
    AllocPlot2 = plt.legend(["Groundwater", "Surface water"], loc='upper left', title="Water source")

    return(AllocPlot1, AllocPlot2)


def allo_plt(df, yaxis_mag=1000000, yaxis_lab='Million', start='2006', end='2015', cat=['tot_allo', 'meter_allo', 'meter_usage'], col_pal='pastel', export_path='', export_name='tot_allo_use_restr.png'):
    """
    Function to plot either total allocation with restrictions or total allo, metered allo, and metered usage with restrictions over a period of years.
    """
    from os import path
    from numpy import in1d
    import seaborn as sns
    from pandas import to_datetime, melt, Period
    import matplotlib.pyplot as plt
#    import matplotlib.patheffects as pathe
#    import matplotlib.ticker as ticker
#    import matplotlib.dates as mdates

    ### Reorganize data
    allo1 = df.sum(axis=0, level=0)
    allo1.index = to_datetime(allo1.index)
    allo1.index.name = 'date'
    allo2 = allo1[start:end] * 1 / yaxis_mag
    dict1 = {'tot_allo': 'tot_allo', 'meter_allo': 'allo', 'meter_usage': 'usage'}
#    dict2 = {'tot_allo': 'tot_ann_up_allo_m3', 'meter_allo': 'ann_up_allo_m3', 'meter_usage': 'usage_m3'}
    lst1 = [dict1[d] for d in cat]
#    lst2 = [dict2[d] for d in cat]

    if allo2.size > 1:
#        index1 = allo2.index.astype('str').str[0:4].astype('int')
#        allo2['index'] = index2

        pw = len(str(yaxis_mag)) - 1

        allo_all = melt(allo2.reset_index(), id_vars='date', value_vars=['tot_allo', 'allo', 'usage'], var_name='tot_allo')
#        allo_up_all = melt(allo2.reset_index(), id_vars='dates', value_vars=['tot_ann_up_allo_m3', 'ann_up_allo_m3', 'usage_m3'], var_name='up_allo')
#        allo_up_all.loc[allo_up_all.up_allo == 'usage_m3', 'value'] = 0

        allo_all = allo_all[in1d(allo_all.tot_allo, lst1)]
#        allo_up_all = allo_up_all[in1d(allo_up_all.up_allo, lst2)]

        index1 = allo_all.date.astype('str').str[0:4].astype('int')
        index2 = [Period(d) for d in index1.tolist()]

        ### Total Allo and restricted allo and usage
        ## Set basic plot settings
        sns.set_style("whitegrid")
        sns.set_context('poster')

        ## Plot total allo
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.barplot(x=index2, y='value', hue='tot_allo', data=allo_all, palette=col_pal, edgecolor='0')
#        sns.barplot(x=index2, y='value', hue='up_allo', data=allo_up_all, palette=col_pal, edgecolor='0', hatch='/')
#        plt.ylabel('Water Volume $(10^{' + str(pw) + '} m^{3}/year$)')
        plt.ylabel('Water Volume $(' + yaxis_lab + '\; m^{3}/year$)')
        plt.xlabel('Year')

        # Legend
        if len(cat) == 3:
            handles, lbs = ax.get_legend_handles_labels()
            order1 = [0,3,1,4,2]
            labels = ['Total allocation', 'Total allocation with restrictions', 'Metered allocation', 'Metered allocation with restrictions', 'Metered usage']
            leg1 = plt.legend([handles[i] for i in order1], labels, loc='upper left')
#        leg1.legendPatch.set_path_effects(pathe.withStroke(linewidth=5, foreground="w"))
        else:
            handles, lbs = ax.get_legend_handles_labels()
            order1 = [0,3]
            plt.legend(handles, ['Total allocation', 'Total allocation with restrictions'], loc='upper left')

        xticks = ax.get_xticks()
        if len(xticks) > 15:
            for label in ax.get_xticklabels()[::2]:
                label.set_visible(False)
            ax.xaxis_date()
            fig.autofmt_xdate(ha='center')
            plt.tight_layout()
        plt.tight_layout()
#      sns.despine(offset=10, trim=True)
        plot2 = ax.get_figure()
        plot2.savefig(path.join(export_path, export_name))

        ### Return plot
        return(ax)


def allo_multi_plot(df, agg_level=[0, 1], index_level=1, plot_fun=allo_plt, export_name='tot_allo_use_restr.png', export_path='', **kwargs):
    """
    Function that calls allo_plt to plot summarized groups from the allo query into several plots.
    """
#    from allo_use_plot_fun import allo_plt

    ### Aggregate to specific level
    allo1 = df.sum(axis=0, level=agg_level)
#    index_name = allo1.index.names[agg_level]
    allo2 = allo1.reset_index()
    index_unique = allo2.ix[:, index_level].unique()

    ### Loop through the index
    for i in index_unique:
        df2 = allo1.xs(i, level=index_level)
        plot_fun(df2, export_name=i + '_' + export_name, export_path=export_path, **kwargs)


def allo_stacked_plt(df, yaxis_mag=1000000, yaxis_lab='Million', start='1990', end='2015', agg_level=[0, 1], cat=['tot_allo'], cat_type='use_type', col_pal='pastel', export_path='', export_name='tot_allo_type.png'):
    """
    Function to plot the summarized groups as stacked bars of the total.
    """
    from os import path
    from numpy import in1d, where, array
    import seaborn as sns
    from pandas import to_datetime, melt, Period
    import matplotlib.pyplot as plt

    ### Reorganize data
    ## Set up dictionaries and parameters
    dict_type = {'public_supply': 'Public Water Supply', 'irrigation': 'Irrigation', 'stockwater': 'Stockwater', 'other': 'Other', 'industry': 'Industry'}
    order1 = array(['public_supply', 'irrigation', 'stockwater', 'other', 'industry'])
    cols1 = sns.color_palette(col_pal)
    col_dict = {'public_supply': cols1[0], 'irrigation': cols1[1], 'stockwater': cols1[2], 'other': cols1[4], 'industry': cols1[3]}

    df2 = df[['tot_ann_allo_m3', 'ann_allo_m3', 'usage_m3']]
    df2 = df2.sum(axis=0, level=agg_level)
    df3 = df2.stack()
    if cat_type is 'use_type':
        df4 = df3.unstack(1)
        cols1 = order1[in1d(order1, df4.columns)]
        df4 = df4[cols1]
        dict0 = dict_type
    else:
        df4 = df3.unstack(1)
    df5 = df4.cumsum(axis=1)
    allo1 = df5.unstack()

    grp1 = df4.columns.tolist()
    lab_names = [dict0[i] for i in grp1]
    col_lab = [col_dict[i] for i in grp1]

    allo1.index = to_datetime(allo1.index)
    allo1.index.name = 'dates'
    allo2 = allo1[start:end] * 1 / yaxis_mag
    dict1 = {'tot_allo': 'tot_ann_allo_m3', 'meter_allo': 'ann_allo_m3', 'meter_usage': 'usage_m3'}
    lst1 = [dict1[d] for d in cat]

#    colors = sns.color_palette(col_pal)

    if allo2.size > 1:
        ### Set plotting parameters
        sns.set_style("whitegrid")
        sns.set_context('poster')
        pw = len(str(yaxis_mag)) - 1

        fig, ax = plt.subplots(figsize=(15, 10))

        for i in grp1[::-1]:
            seq1 = int(where(in1d(grp1, i))[0])
            allo_all = melt(allo2[i].reset_index(), id_vars='dates', value_vars=['tot_ann_allo_m3', 'ann_allo_m3', 'usage_m3'], var_name=i)
            allo_all = allo_all[in1d(allo_all[i], lst1)]

            index1 = allo_all.dates.astype('str').str[0:4].astype('int')
            index2 = [Period(d) for d in index1.tolist()]
            sns.barplot(x=index2, y='value', data=allo_all, edgecolor='0', color=col_lab[seq1], label=i)

#        plt.ylabel('Allocated Water Volume $(10^{' + str(pw) + '} m^{3}/year$)')
        plt.ylabel('Water Volume $(' + yaxis_lab + '\; m^{3}/year$)')
        plt.xlabel('Year')

        # Legend
        handles, lbs = ax.get_legend_handles_labels()
        plt.legend(handles, lab_names[::-1], loc='upper left')

        xticks = ax.get_xticks()
        if len(xticks) > 15:
            for label in ax.get_xticklabels()[::2]:
                label.set_visible(False)
            ax.xaxis_date()
            fig.autofmt_xdate(ha='center')
            plt.tight_layout()
        plt.tight_layout()
#      sns.despine(offset=10, trim=True)
        plot2 = ax.get_figure()
        plot2.savefig(path.join(export_path, export_name))

        ### Return plot
        return(ax)


def allo_restr_plt(df, yaxis_mag=1000000, yaxis_lab='Million', start='2006', end='2015', cat=['tot_allo', 'meter_allo', 'meter_usage'], col_pal='pastel', export_path='', export_name='tot_allo_use_restr.png'):
    """
    Function to plot either total allocation with restrictions or total allo, metered allo, and metered usage with restrictions over a period of years.
    """
    from os import path
    from numpy import in1d
    import seaborn as sns
    from pandas import to_datetime, melt, Period
    import matplotlib.pyplot as plt
    from collections import OrderedDict
#    import matplotlib.patheffects as pathe
#    import matplotlib.ticker as ticker
#    import matplotlib.dates as mdates

    base_cat = ['tot_allo', 'meter_allo', 'meter_usage']

    ### Reorganize data
    allo1 = df.sum(axis=0, level=0)
    allo1.index = to_datetime(allo1.index)
    allo1.index.name = 'date'
    allo2 = allo1[start:end] * 1 / yaxis_mag
    dict1 = {'tot_allo': 'tot_allo', 'meter_allo': 'allo', 'meter_usage': 'usage'}
    dict2 = {'tot_allo': 'tot_allo_restr', 'meter_allo': 'allo_restr', 'meter_usage': 'usage'}
    lst1 = [dict1[d] for d in cat]
    lst2 = [dict2[d] for d in cat]

    if allo2.size > 1:
#        index1 = allo2.index.astype('str').str[0:4].astype('int')
#        allo2['index'] = index2

        pw = len(str(yaxis_mag)) - 1

        allo_all = melt(allo2.reset_index(), id_vars='date', value_vars=['tot_allo', 'allo', 'usage'], var_name='tot_allo')
        allo_up_all = melt(allo2.reset_index(), id_vars='dates', value_vars=['tot_allo_restr', 'allo_restr', 'usage'], var_name='up_allo')
        allo_up_all.loc[allo_up_all.up_allo == 'usage', 'value'] = 0

        allo_all = allo_all[in1d(allo_all.tot_allo, lst1)]
        allo_up_all = allo_up_all[in1d(allo_up_all.up_allo, lst2)]

        index1 = allo_all.date.astype('str').str[0:4].astype('int')
        index2 = [Period(d) for d in index1.tolist()]

        ### Total Allo and restricted allo and usage
        ## Set basic plot settings
        sns.set_style("whitegrid")
        sns.set_context('poster')
        col_pal1 = sns.color_palette(col_pal)
        col_pal2 = [col_pal1[base_cat.index(i)] for i in cat]

        ## Plot total allo
        fig, ax = plt.subplots(figsize=(15, 10))
        sns.barplot(x=index2, y='value', hue='tot_allo', data=allo_all, palette=col_pal2, edgecolor='0')
        sns.barplot(x=index2, y='value', hue='up_allo', data=allo_up_all, palette=col_pal2, edgecolor='0', hatch='/')
#        plt.ylabel('Water Volume $(10^{' + str(pw) + '} m^{3}/year$)')
        plt.ylabel('Water Volume $(' + yaxis_lab + '\; m^{3}/year$)')
        plt.xlabel('Year')

        # Legend
        handles, lbs = ax.get_legend_handles_labels()
        hand_len = len(handles)
        order1 = [lbs.index(j) for j in ['tot_allo', 'tot_allo_restr', 'allo', 'allo_restr', 'usage'] if j in lbs]
        label_dict1 = OrderedDict([('tot_allo', ['Total allocation', 'Total allocation with restrictions']), ('meter_allo', ['Metered allocation', 'Metered allocation with restrictions']), ('meter_usage', ['Metered usage'])])
        labels = []
        [labels.extend(label_dict1[i]) for i in cat]
        leg1 = plt.legend([handles[i] for i in order1], labels, loc='upper left')
#        leg1.legendPatch.set_path_effects(pathe.withStroke(linewidth=5, foreground="w"))

        xticks = ax.get_xticks()
        if len(xticks) > 15:
            for label in ax.get_xticklabels()[::2]:
                label.set_visible(False)
            ax.xaxis_date()
            fig.autofmt_xdate(ha='center')
            plt.tight_layout()
        plt.tight_layout()
#      sns.despine(offset=10, trim=True)
        plot2 = ax.get_figure()
        plot2.savefig(path.join(export_path, export_name))

        ### Return plot
        return(ax)













