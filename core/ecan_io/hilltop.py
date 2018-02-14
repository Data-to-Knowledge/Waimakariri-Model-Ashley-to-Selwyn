# -*- coding: utf-8 -*-
"""
Hilltop read functions.
Hilltop uses a fixed base date as 1940-01-01, while the standard unix/POSIT base date is 1970-01-01.
"""
# import Hilltop
import os
from configparser import ConfigParser
from win32com.client import Dispatch, pywintypes
# from lxml import etree
import pandas as pd
import numpy as np
import geopandas as gpd
from core.ecan_io.mssql import rd_sql
from core.misc import select_sites, pytime_to_datetime, time_switch
from core.spatial.vector import xy_to_gpd, sel_sites_poly
from datetime import datetime

#####################################################
#### New method - not ready yet...


# def ht_sites(hts, sites=None):
#     """
#     Function to read all of the sites in an hts file and the associated site info.
#
#     hts -- Path to hts file (str).\n
#     sites -- Optional list, array, series of site names to return.
#     """
#
#     dfile1 = Hilltop.Connect(hts)
#     site_lst = Hilltop.SiteList(dfile1)
#
#     if not site_lst:
#         print('No sites in ' + hts)
#         return(DataFrame())
#
#     if sites is not None:
#         sites1 = select_sites(sites)
#         site_lst = [i for i in site_lst if i in sites1]
#
#     site_info = DataFrame()
#
#     for i in site_lst:
#         try:
#             info1 = Hilltop.MeasurementList(dfile1, i)
#         except SystemError:
# #            print('Site ' + str(i) + " didn't work")
#             continue
#         info1.loc[:, 'site'] = i
#         site_info = concat([site_info, info1])
#     site_info.reset_index(drop=True, inplace=True)
#
#     site_info.loc[:, 'Start Time'] = to_datetime(site_info.loc[:, 'Start Time'], format='%d-%b-%Y %H:%M:%S')
#     site_info.loc[:, 'End Time'] = to_datetime(site_info.loc[:, 'End Time'], format='%d-%b-%Y %H:%M:%S')
#
#     len_all = len(site_lst)
#     len_got = len(site_info.site.unique())
#     print('Missing ' + str(len_all - len_got) + ' sites, which is ' + str(round(100 * ((len_all - len_got)/len_all), 1)) + '% of the total')
#
#     Hilltop.Disconnect(dfile1)
#     return(site_info)
#
#
# def ht_get_data(hts, sites=None, from_date=None, to_date=None, agg_method='Average', interval='1 day', alignment='00:00', output_missing_sites=False):
#     """
#
#     """
#
#     site_info = ht_sites(hts, sites=sites)
#
#     if isinstance(from_date, str):
#         from_date1 = to_datetime(from_date).strftime('%d-%b-%Y %H:%M')
#     else:
#         from_date1 = ''
#     if isinstance(to_date, str):
#         to_date1 = to_datetime(to_date).strftime('%d-%b-%Y %H:%M')
#     else:
#         to_date1 = ''
#
#     data1 = []
#     missing_site_data = DataFrame()
#
#     dfile1 = Hilltop.Connect(hts)
#     for i in site_info.index:
#         site = site_info.loc[i, 'site']
#         mtype = site_info.loc[i, 'Measurement']
#
#         d1 = Hilltop.GetData(dfile1, site, mtype, from_date1, to_date1, method=agg_method, interval=interval, alignment=alignment)
#
#         if d1.empty:
#             missing_site_data = concat([missing_site_data, site_info.loc[i]])
# #            print('No data for site ' + str(site) + ' and mtype ' + str(mtype))
#             continue
#
#         if (interval == '1 day') & (alignment == '00:00'):
#             d1.index = d1.index.normalize()
#         d1.name = 'data'
#         d1.index.name = 'time'
#         d2 = d1.reset_index()
#         d2.loc[:, 'site'] = site
#         d2.loc[:, 'mtype'] = mtype
#         print(site, mtype)
#
#         data1.append(d2)
#
#     try:
#         data2 = concat(data1)
#     except MemoryError:
#         print('Not enough memory for a 32bit application')
#
#     if missing_site_data.empty:
#         print('No missing data for any sites/mtypes')
#     else:
#         print('No data for ' + str(len(missing_site_data)) + ' sites/mtype combos')
#         if output_missing_sites:
#             return(data2, missing_site_data)
#     return(data2)


######################################################
#### Old access method


def parse_dsn(dsn_path):
    """
    Function to parse a dsn file and all sub-dsn files into paths to hts files. Returns a list of hts paths.
    """

    base_path = os.path.dirname(dsn_path)

    dsn = ConfigParser()
    dsn.read(dsn_path)
    files1 = [os.path.join(base_path, i[1]) for i in dsn.items('Hilltop') if 'file' in i[0]]
    hts1 = [i for i in files1 if i.endswith('.hts')]
    dsn1 = [i for i in files1 if i.endswith('.dsn')]
    while dsn1:
        for f in dsn1:
            base_path = os.path.dirname(f)
            p1 = ConfigParser()
            p1.read(f)
            files1 = [os.path.join(base_path, i[1]) for i in p1.items('Hilltop') if 'file' in i[0]]
            hts1.extend([i for i in files1 if i.endswith('.hts')])
            dsn1.remove(f)
            dsn1[0:0] = [i for i in files1 if i.endswith('.dsn')]
    return(hts1)


def rd_hilltop_sites(hts, sites=None, mtypes=None, rem_wq_sample=True):
    """
    Function to read the site names, measurement types, and units of a Hilltop hts file. Returns a DataFrame.

    Arguments:\n
    hts -- Path to the hts file (str).\n
    sites -- A list of site names within the hts file.\n
    mtypes -- A list of measurement types that should be returned.
    """

    if sites is not None:
        sites = select_sites(sites)
    if mtypes is not None:
        mtypes = select_sites(mtypes)

    cat = Dispatch("Hilltop.Catalogue")
    if not cat.Open(hts):
        raise ValueError(cat.errmsg)

    dfile = Dispatch("Hilltop.DataRetrieval")
    try:
        dfile.Open(hts)
    except ValueError:
        print(dfile.errmsg)

    sites_lst = []

    ### Iterate through all sites/datasources/mtypes
    cat.StartSiteEnum
    while cat.GetNextSite:
        site_name = cat.SiteName
        if sites is None:
            pass
        elif site_name in sites:
            pass
        else:
            continue
        while cat.GetNextDataSource:
            ds_name = cat.DataSource
            try:
                start1 = pytime_to_datetime(cat.DataStartTime)
                end1 = pytime_to_datetime(cat.DataEndTime)
            except ValueError:
                bool_site = dfile.FromSite(site_name, ds_name, 1)
                if bool_site:
                    start1 = pytime_to_datetime(cat.DataStartTime)
                    end1 = pytime_to_datetime(cat.DataEndTime)
                else:
                    print('No site data for ' + site_name + '...for some reason...')
            while cat.GetNextMeasurement:
                mtype1 = cat.Measurement
                if mtype1 == 'Item2':
                    continue
                elif mtypes is None:
                    pass
                elif mtype1 in mtypes:
                    pass
                else:
                    continue
                divisor = cat.Divisor
                unit1 = cat.Units
                if unit1 == '%':
#                    print('Site ' + name1 + ' has no units')
                    unit1 = ''
                sites_lst.append([site_name, ds_name, mtype1, unit1, divisor, str(start1), str(end1)])

    sites_df = pd.DataFrame(sites_lst, columns=['site', 'data_source', 'mtype', 'unit', 'divisor', 'start_date', 'end_date'])
    if rem_wq_sample:
        sites_df = sites_df[~(sites_df.mtype == 'WQ Sample')]
    dfile.Close()
    cat.Close()
    return(sites_df)


def rd_ht_quan_data(hts, sites=None, mtypes=None, start=None, end=None, agg_period=None, agg_n=1, fun=None, output_site_data=False, exclude_mtype=None):
    """
    Function to read data from an hts file and optionally select specific sites and aggregate the data.

    Parameters
    ----------
    hts : str
        Path to the hts file.
    sites : list
        A list of site names within the hts file.
    mtypes : list
        A list of measurement types that should be returned.
    start : str
        The start date to retreive from the data in ISO format (e.g. '2011-11-30 00:00').
    end : str
        The end date to retreive from the data in ISO format (e.g. '2011-11-30 00:00').
    agg_period : str
        The resample period (e.g. 'day', 'month').
    agg_n : int
        The number of periods (e.g. 1 for 1 day).
    fun : str
        The resampling function.
    output_site_data : bool
        Should the sites data be output?

    Returns
    -------
    DataFrame
    """

    agg_name_dict = {'sum': 4, 'count': 5, 'mean': 1}
    agg_unit_dict = {'l/s': 1, 'm3/s': 1, 'm3/hour': 1, 'mm': 1, 'm3': 4}
    unit_convert = {'l/s': 0.001, 'm3/s': 1, 'm3/hour': 1, 'mm': 1, 'm3': 1}

    ### First read all of the sites in the hts file and select the ones to be read
    sites_df = rd_hilltop_sites(hts, sites=sites, mtypes=mtypes)
    sites_df = sites_df[sites_df.unit.isin(list(agg_unit_dict.keys()))]
    if isinstance(exclude_mtype, list):
        sites_df = sites_df[~sites_df.mtype.isin(exclude_mtype)]

    ### Select out the sites/mtypes within the date range
    if isinstance(start, str):
        sites_df = sites_df[sites_df.end_date >= start]
    if isinstance(end, str):
        sites_df = sites_df[sites_df.start_date <= end]

    ### Open the hts file
    dfile = Dispatch("Hilltop.DataRetrieval")
    try:
        dfile.Open(hts)
    except ValueError:
        print(dfile.errmsg)

    ### Iterate through he hts file
    df_lst = []
    for i in sites_df.index:
        site = sites_df.loc[i, 'site']
#        data_source = sites_df.loc[i, 'data_source']
        mtype = sites_df.loc[i, 'mtype']
        unit = sites_df.loc[i, 'unit']
#        if mtype == 'Flow':
#            mtype = 'Flow [Flow]'
        if fun is None:
            agg_val = agg_unit_dict[unit]
        else:
            agg_val = agg_name_dict[fun]
        if dfile.FromSite(site, mtype, 1):

            ## Set up start and end times and aggregation initiation
            start_time = pytime_to_datetime(dfile.DataStartTime)
            end_time = pytime_to_datetime(dfile.DataEndTime)
            if (start_time.year < 1900) | (end_time.year < 1900):
                print('Site ' + site + ' has a start or end time prior to 1900')
                continue
            if (start is None):
                if (agg_period is not None):
                    start1 = str(pd.to_datetime(start_time).ceil(str(agg_n) + time_switch(agg_period)).date())
                else:
                    start1 = dfile.DataStartTime
            else:
                start1 = start
            if end is None:
                end1 = dfile.DataEndTime
            else:
                end1 = end
            if not dfile.FromTimeRange(start1, end1):
                continue
            if (agg_period is not None):
                dfile.SetMode(agg_val, str(agg_n) + ' ' + agg_period)

            ## Extract data
            data = []
            time = []
            if dfile.getsinglevbs == 0:
                t1 = dfile.value
                if isinstance(t1, str):
                    print('site ' + site + ' has nonsense data')
                else:
                    data.append(t1)
                    time.append(str(pytime_to_datetime(dfile.time)))
                    while dfile.getsinglevbs == 0:
                        data.append(dfile.value)
                        time.append(str(pytime_to_datetime(dfile.time)))
                    if data:
                        df_temp = pd.DataFrame({'time': time, 'data': data, 'site': site, 'mtype': mtype})
                        df_lst.append(df_temp)

    dfile.Close()
    if df_lst:
        df1 = pd.concat(df_lst)
        df1.loc[:, 'time'] = pd.to_datetime(df1.loc[:, 'time'])
        df2 = df1.set_index(['mtype', 'site', 'time']).data * unit_convert[unit]
    else:
        df2 = pd.DataFrame([], index=['mtype', 'site', 'time'])
    if output_site_data:
        return df2, sites_df
    else:
        return df2


def rd_ht_wq_data(hts, sites=None, mtypes=None, start=None, end=None, dtl_method=None, output_site_data=False, mtype_params=None, sample_params=None):
    """
    Function to read data from an hts file and optionally select specific sites and aggregate the data.

    Parameters
    ----------
    hts : str
        Path to the hts file.
    sites : list
        A list of site names within the hts file.
    mtypes : list
        A list of measurement types that should be returned.
    start : str
        The start date to retreive from the data in ISO format (e.g. '2011-11-30 00:00').
    end : str
        The end date to retreive from the data in ISO format (e.g. '2011-11-30 00:00').
    dtl_method : None, 'standard', 'trend'
        The method to use to convert values under a detection limit to numeric. None does no conversion. 'standard' takes half of the detection limit. 'trend' is meant as an output for trend analysis with includes an additional column dtl_ratio referring to the ratio of values under the detection limit.
    output_site_data : bool
        Should the site data be output?

    Returns
    -------
    DataFrame
    """

#    agg_unit_dict = {'l/s': 1, 'm3/s': 1, 'm3/hour': 1, 'mm': 1, 'm3': 4}
#    unit_convert = {'l/s': 0.001, 'm3/s': 1, 'm3/hour': 1, 'mm': 1, 'm3': 4}

    sites1 = select_sites(sites)

    #### Extract by polygon
    if isinstance(sites1, gpd.GeoDataFrame):
        ## Surface water sites
        sw_sites_tab = rd_sql('SQL2012PROD05', 'Squalarc', 'SITES', col_names=['SITE_ID', 'NZTMX', 'NZTMY'])
        sw_sites_tab.columns = ['site', 'NZTMX', 'NZTMY']
        gdf_sw_sites = xy_to_gpd('site', 'NZTMX', 'NZTMY', sw_sites_tab)
        sites1a = sites1.to_crs(gdf_sw_sites.crs)
        sw_sites2 = sel_sites_poly(gdf_sw_sites, sites1a).drop('geometry', axis=1)

        ## Groundwater sites
        gw_sites_tab = rd_sql('SQL2012PROD05', 'Wells', 'WELL_DETAILS', col_names=['WELL_NO', 'NZTMX', 'NZTMY'])
        gw_sites_tab.columns = ['site', 'NZTMX', 'NZTMY']
        gdf_gw_sites = xy_to_gpd('site', 'NZTMX', 'NZTMY', gw_sites_tab)
        gw_sites2 = sel_sites_poly(gdf_gw_sites, sites1a).drop('geometry', axis=1)

        sites2 = sw_sites2.site.append(gw_sites2.site).astype(str).tolist()
    else:
        sites2 = sites1

    ### First read all of the sites in the hts file and select the ones to be read
    sites_df = rd_hilltop_sites(hts, sites=sites2, mtypes=mtypes, rem_wq_sample=False)

    ### Select out the sites/mtypes within the date range
    if isinstance(start, str):
        sites_df = sites_df[sites_df.end_date >= start]
    if isinstance(end, str):
        sites_df = sites_df[sites_df.start_date <= end]

    ### Open the hts file
    wqr = Dispatch("Hilltop.WQRetrieval")
    dfile = Dispatch("Hilltop.DataFile")
    try:
        dfile.Open(hts)
    except ValueError:
        print(dfile.errmsg)

    ### Iterate through he hts file
    df_lst = []
    for i in sites_df.index:
        site_data = sites_df.loc[i]
        site = site_data['site']
        mtype = site_data['mtype']
        if mtype == 'WQ Sample':
            continue
        wqr = dfile.FromWQSite(site, mtype)

        ## Set up start and end times and aggregation initiation
        if start is None:
            start1 = wqr.DataStartTime
        else:
            start1 = pywintypes.TimeType.strptime(start, '%Y-%m-%d')
        if end is None:
            end1 = wqr.DataEndTime
        else:
            end1 = pywintypes.TimeType.strptime(end, '%Y-%m-%d')

        if not wqr.FromTimeRange(start1, end1):
            continue

        ## Extract data
        data = []
        time = []
        sample_p = []

        test_params = sites_df[sites_df.site == site].mtype.unique()
        if ('WQ Sample' in test_params) & (isinstance(mtype_params, list) | isinstance(sample_params, list)):
            mtype_p = []
            while wqr.GetNext:
                data.append(wqr.value)
                time.append(str(pytime_to_datetime(wqr.time)))
                sample_p.append({sp: wqr.params(sp).encode('ascii', 'ignore').decode() for sp in sample_params})
                mtype_p.append({mp: wqr.params(mp).encode('ascii', 'ignore').decode() for mp in mtype_params})
        else:
            while wqr.GetNext:
                data.append(wqr.value)
                time.append(str(pytime_to_datetime(wqr.time)))

        if data:
            df_temp = pd.DataFrame({'time': time, 'data': data, 'site': site, 'mtype': mtype})
            if sample_p:
                df_temp = pd.concat([df_temp, pd.DataFrame(sample_p), pd.DataFrame(mtype_p)], axis=1)
            df_lst.append(df_temp)

    dfile.Close()
    wqr.close()
    if df_lst:
        data = pd.concat(df_lst)
        data.loc[:, 'time'] = pd.to_datetime(data.loc[:, 'time'])
        data1 = pd.to_numeric(data.loc[:, 'data'], errors='coerce')
        data.loc[data1.notnull(), 'data'] = data1[data1.notnull()]
#        data.loc[:, 'data'].str.replace('*', '')
        data = data.reset_index(drop=True)

        #### Convert detection limit values
        if dtl_method is not None:
            less1 = data['data'].str.match('<')
            if less1.sum() > 0:
                less1.loc[less1.isnull()] = False
                data2 = data.copy()
                data2.loc[less1, 'data'] = pd.to_numeric(data.loc[less1, 'data'].str.replace('<', ''), errors='coerce') * 0.5
                if dtl_method == 'standard':
                    data3 = data2
                if dtl_method == 'trend':
                    df1 = data2.loc[less1]
                    count1 = data.groupby('mtype')['data'].count()
                    count1.name = 'tot_count'
                    count_dtl = df1.groupby('mtype')['data'].count()
                    count_dtl.name = 'dtl_count'
                    count_dtl_val = df1.groupby('mtype')['data'].nunique()
                    count_dtl_val.name = 'dtl_val_count'
                    combo1 = pd.concat([count1, count_dtl, count_dtl_val], axis=1, join='inner')
                    combo1['dtl_ratio'] = (combo1['dtl_count'] / combo1['tot_count']).round(2)

                    ## conditionals
                    param2 = combo1[(combo1['dtl_ratio'] > 0.4) & (combo1['dtl_val_count'] != 1)]
                    over_40 = data['mtype'].isin(param2.index)

                    ## Calc detection limit values
                    data3 = pd.merge(data, combo1['dtl_ratio'].reset_index(), on='mtype', how='left')
                    data3.loc[:, 'data_dtl'] = data2['data']

                    max_dtl_val = data2[over_40 & less1].groupby('mtype')['data'].transform('max')
                    max_dtl_val.name = 'dtl_data_max'
                    data3.loc[over_40 & less1, 'data_dtl'] = max_dtl_val
            else:
                data3 = data
        else:
            data3 = data

        if output_site_data:
            sites_df = sites_df[~(sites_df.mtype == 'WQ Sample')]
            return data3, sites_df
        else:
            return data3


# def rd_ht_xml_sites(xml):
#     """
#     Function to read a Hilltop xml file and return the site names. The xml file should be a complete export of an hts file.
#     """
#
#     ### Parse xml
#     root = etree.iterparse(xml, tag='Measurement')
#
#     ### Iterate
#     sites = []
#     mtypes = []
#     for event, elem in root:
#         sites.append(elem.values()[0])
#         ds = elem.find('DataSource')
#         mtypes.append(ds.values()[0])
#         elem.clear()
#         while elem.getprevious() is not None:
#             del elem.getparent()[0]
#
#     ### Return
#     df = pd.DataFrame([sites, mtypes]).transpose()
#     df.columns = ['site', 'mtype']
#     return(df)
#
#
# def parse_ht_xml(xml, ht_fun, select=None, corr_csv=r'C:\ecan\hilltop\ht_corrections.csv'):
#     """
#     Function to read a Hilltop xml file and apply a function on each individual site time series. The input to the function is a single pandas time series. The output should be a Series or DataFrame. Specific sites with specific mtypes can be passed in the form of a two column DataFrame with headers as 'site' and 'mtype'.
#     """
#
#     ### Base parameters
#     rem_s = 10958*24*60*60
#     corr = pd.read_csv(corr_csv)
#     xml_name = basename(xml)
#
#     ### Select corrections
#     corr1 = corr[corr.file_name == xml_name]
#
#     ### Parse xml
#     root = etree.iterparse(xml, tag='Measurement')
#
#     ### Iterate
#     results1 = []
#     for event, elem in root:
#         ## Get data
#         site = elem.values()[0]
#         mtype = elem.find('DataSource').values()[0]
#
#         if (select is not None):
#             if (not isinstance(select, DataFrame)):
#                 raise ValueError('Make sure the input is a DataFrame with two columns!')
#             elif all(select.columns == ['site', 'mtype']):
#                 site_check = any([set([site, mtype]) == set([select.loc[i].site, select.loc[i].mtype]) for i in select.index])
#                 if not site_check:
#                     continue
#
# #        units = elem.find('DataSource').find('ItemInfo').find('Units').text
#         site_data = [j.text.split() for j in elem.find('Data').findall('V')]
#
#         ## Convert to dataframe
#         o2 = DataFrame(site_data, columns=['date', 'val'])
#         o2.loc[:,['date', 'val']] = o2.loc[:,['date', 'val']].astype(float)
#         o2.loc[:, 'date'] = to_datetime(o2.loc[:, 'date'] - rem_s, unit='s')
#         o2.set_index('date', inplace=True)
#
#         ## Make corrections
#         corr_index = (corr1.orig_site == site) & (corr1.orig_mtype == mtype)
#         if any(corr_index):
#             site, mtype = corr1.loc[corr_index, ['new_site', 'new_mtype']].values.tolist()[0]
#
#         ## Clear element from memory
#         elem.clear()
#         while elem.getprevious() is not None:
#             del elem.getparent()[0]
#
#         ## Do stats
#         stats1 = ht_fun(o2, mtype, site)
#
#         ## Add additional site specific columns/data
# #        stats1.loc[:, 'site'] = site
# #        stats1.loc[:, 'mtype'] = mtype
# #        stats1.loc[:, 'units'] = units
#
#         ## Append
#         results1.append(stats1)
#
#     ### Combine data
#     df_out = concat(results1)
#
#     ### Return
#     return(df_out)


# def data_check_fun(data, mtype, site):
#     """
#     Various data checks on the hilltop data. This function should be an input to parse_ht_xml.
#     """
#
#     def infer_freq1(x):
#         if len(x) > 7:
#             freq1 = infer_freq(x.index[3:6], warn=False)
#         else:
#             freq1 = None
#         return(freq1)
#
#     ## Run stat checks
#     first1 = data.index.min()
#     last1 = data.index.max()
#     n_neg = (data < 0).sum()[0]
#     n_zero = (data == 0).sum()[0]
#     count1 = data.count()[0]
#     freq = infer_freq1(data)
#     neg_ratio = round(n_neg/float(count1), 2)
#     zero_ratio = round(n_zero/float(count1), 2)
#
#     ## Construct output
#     out1 = [first1, last1, freq, count1, n_neg, n_zero, neg_ratio, zero_ratio]
#     out_names = ['start_date', 'end_date', 'time_res', 'n_data', 'n_neg', 'n_zero', 'n_neg/n_data', 'n_zero/n_data']
#     df_out = DataFrame([out1], columns=out_names)
#     df_out.loc[:, 'site'] = site
#     df_out.loc[:, 'mtype'] = mtype
#
#     ### return
#     return(df_out)
#
#
# def iter_xml_dir(fpath, stats_fun, with_xml=False, select=None, export=False, export_name='results.csv'):
#
#     ### Read files in directory
#     files = rd_dir(fpath, 'xml')
#
#     ### Iterate through each file
#     list_out = []
#     for i in files:
#         print(i)
#         if select is None:
#             out1 = parse_ht_xml(join(fpath, i), stats_fun)
#         elif isinstance(select, DataFrame) & (len(select.columns) == 3):
#             select1 = select.loc[select.file_name == i, ['site', 'mtype']]
#             out1 = parse_ht_xml(join(fpath, i), stats_fun, select1)
#         if with_xml:
#             out1['xml'] = [i] * len(out1)
#         list_out.append(out1)
#
#     ### Combine
#     df_out = concat(list_out)
#     df_out.index.name = 'date'
#
#     ### Export and return
#     if export:
#         df_out.to_csv(join(fpath, export_name), encoding='utf-8')
#     return(df_out)
#
#
# def all_data_fun(data, mtype, site):
#     """
#     Function for parse_ht_xml to return all data.
#     """
#     data.loc[:, 'site'] = site
#     data.loc[:, 'mtype'] = mtype
#     return(data)
#
#
# def proc_use_data(data, mtype, site, time_period='D', n_std=4):
#     """
#     Function for parse_ht_xml to process the data and aggregate it to a defined resolution.
#     """
#
#     ### Select the process sequence based on the mtype and convert to period volume
#     data[data < 0] = nan
#     count1 = float(data.count().values[0])
#
#     if mtype == 'Water Meter':
#         ## Check to determine whether it is cumulative or period volume
#         diff1 = data.diff()[1:]
#         neg_index = diff1 < 0
#         neg_ratio = sum(neg_index.values)/count1
#         if neg_ratio > 0.1:
#             outliers = abs(data - data.mean()) > (data.std() * n_std)
#             data[outliers] = nan
#             vol = data
#         else:
#             # Replace the negative values with zero and the very large values
#             diff1[diff1 < 0] = data[diff1 < 0]
#             outliers = abs(diff1 - diff1.mean()) > (diff1.std() * n_std)
#             diff1[outliers] = nan
#             vol = diff1
#     elif (mtype == 'Abstraction Volume') | (mtype == 'Average Flow'):
#         outliers = abs(data - data.mean()) > (data.std() * n_std)
#         data[outliers] = nan
#         vol = data
#     elif mtype == 'Flow':
#         outliers = abs(data - data.mean()) > (data.std() * n_std)
#         data[outliers] = nan
#
#         # Determine the diff index
#         t1 = Series(data.index).diff().dt.seconds.shift(-1)
#         t1.iloc[-1] = t1.iloc[-2]
#         t1.index = data.index
#         # Convert to volume
#         vol = data.multiply(t1, axis=0) * 0.001
#
#     ## Estimate the NAs
#     vol2 = vol.fillna(method='ffill')
#
#     ## Resample the volumes
#     vol_res = vol2.resample(time_period).sum()
#     vol_res.loc[:, 'site'] = site
#
#     return(vol_res)


def convert_site_names(names, rem_m=True):

    names1 = names.str.replace('[:\.]', '/')
#    names1.loc[names1 == 'L35183/580-M1'] = 'L35/183/580-M1' What to do with this one?
    names1.loc[names1 == 'L370557-M1'] = 'L37/0557-M1'
    names1.loc[names1 == 'L370557-M72'] = 'L37/0557-M72'
    names1.loc[names1 == 'BENNETT K38/0190-M1'] = 'K38/0190-M1'
    names1 = names1.str.upper()
    if rem_m:
        list_names1 = names1.str.findall('[A-Z]+\d\d/\d\d\d\d')
        names_len_bool = list_names1.apply(lambda x: len(x)) == 1
        names2 = names1.copy()
        names2[names_len_bool] = list_names1[names_len_bool].apply(lambda x: x[0])
        names2[~names_len_bool] = np.nan
    else:
        list_names1 = names1.str.findall('[A-Z]+\d\d/\d\d\d\d\s*-\s*M\d*')
        names_len_bool = list_names1.apply(lambda x: len(x)) == 1
        names2 = names1.copy()
        names2[names_len_bool] = list_names1[names_len_bool].apply(lambda x: x[0])
        names2[~names_len_bool] = np.nan

    return names2


def proc_ht_use_data(ht_data, n_std=4):
    """
    Function to process the water usage data at daily resolution.
    """

    ### Groupby mtypes and sites
    grp = ht_data.groupby(level=['mtype', 'site'])

    res1 = []
    for index, data1 in grp:
        data = data1.copy()
        mtype, site = index
#        units = sites[(sites.site == site) & (sites.mtype == mtype)].unit

        ### Select the process sequence based on the mtype and convert to period volume

        data[data < 0] = np.nan

        if mtype == 'Water Meter':
            ## Check to determine whether it is cumulative or period volume
            count1 = float(data.count())
            diff1 = data.diff()
            neg_index = diff1 < 0
            neg_ratio = sum(neg_index.values)/count1
            if neg_ratio > 0.1:
                outliers = abs(data - data.mean()) > (data.std() * n_std)
                data[outliers] = np.nan
                vol = data
            else:
                # Replace the negative values with zero and the very large values
                diff1[diff1 < 0] = data[diff1 < 0]
                outliers = abs(diff1 - diff1.mean()) > (diff1.std() * n_std)
                diff1.loc[outliers] = np.nan
                vol = diff1
        elif mtype in ['Compliance Volume', 'Volume']:
            outliers = abs(data - data.mean()) > (data.std() * n_std)
            data.loc[outliers] = np.nan
            vol = data
        elif mtype == 'Flow':
            outliers = abs(data - data.mean()) > (data.std() * n_std)
            data.loc[outliers] = np.nan

#            # Determine the diff index
#            t1 = Series(data.index).diff().dt.seconds.shift(-1)
#            t1.iloc[-1] = t1.iloc[-2]
#            t1.index = data.index
#            # Convert to volume
#            vol = data.multiply(t1, axis=0) * 0.001
            vol = (data * 60*60*24).fillna(method='ffill').round(4)
        elif mtype == 'Average Flow':
            outliers = abs(data - data.mean()) > (data.std() * n_std)
            data.loc[outliers] = np.nan
            vol = (data * 24).fillna(method='ffill').round(4)

        res1.append(vol)

    ### Convert to dataframe
    df1 = pd.concat(res1).reset_index()

    ### Drop the mtypes level and uppercase the sites
    df2 = df1.drop('mtype', axis=1)
    df2.loc[:, 'site'] = df2.loc[:, 'site'].str.upper()

    ### Remove duplicate WAPs
    df3 = df2.groupby(['site', 'time']).data.last()

    return df3

