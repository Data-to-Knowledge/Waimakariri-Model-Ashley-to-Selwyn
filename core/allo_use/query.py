# -*- coding: utf-8 -*-
"""
Functions to query the allocation and usage data.
"""


def allo_query(shp=None, grp_by=['date', 'take_type', 'use_type'], allo_col=['allo'], use_col=['usage'], agg_yr=True, crc='all', crc_rem='none', wap_rem='none', wap='all', take_type='all', use_type='all', catch_num='all', swaz='all', cwms_zone='all', swaz_grp='all', years='all', gr_than='all', gwaz='all', catch_name='all', catch_grp_name='all', from_date=None, to_date=None, sd_only=False, allo_use_file='S:/Surface Water/shared/base_data/usage/allo_est_use_mon.h5', allo_gis_file=r'S:\Surface Water\shared\GIS_base\vector\allocations\allo_gis.shp', allo_restr_file='S:/Surface Water/shared/base_data/usage/allo_use_ros_mon.csv', export=True, export_path='summary1.csv', debug=False):
    """
    Function to query the water use and allocation results data. Allows for the selection/filtering and aggregation of many imbedded fields. Create a list only when you want to filter the data. Otherwise, leave the srguments default.
    """
    from pandas import concat, read_csv, read_hdf, merge, to_numeric, to_datetime
    from numpy import floor
    from core.ts.ts import grp_ts_agg
    from geopandas import read_file
    from core.spatial.vector import sel_sites_poly

    ### Read in the data
#    data = read_hdf(allo_use_file)[['crc', 'dates', 'take_type', 'use_type', 'mon_vol', 'up_allo_m3', 'usage', 'usage_est']]
    data = read_hdf(allo_use_file).drop(['ann_allo_m3', 'ann_usage_m3', 'band', 'band_restr', 'gauge_num', 'ann_restr_allo_m3', 'usage_ratio_est'], axis=1)
    data.rename(columns={'mon_restr_allo_m3': 'allo_restr', 'mon_usage_m3': 'usage', 'mon_allo_m3': 'allo'}, inplace=True)
    allo_gis = read_file(allo_gis_file)

    if shp is not None:
        poly = read_file(shp)
        allo_gis = sel_sites_poly(poly, allo_gis)

    ### Aggregate if needed
    if agg_yr:
        data = grp_ts_agg(data, ['crc', 'take_type', 'allo_block', 'wap'], 'date', 'A-JUN').sum().reset_index()
    df = merge(data, allo_gis.drop(['geometry'] , axis=1), on=['crc', 'take_type', 'allo_block', 'wap'])
    df.loc[:, 'catch_grp'] = to_numeric(df.loc[:, 'catch_grp'], errors='coerse')

    ### Run through filters
    if type(crc) == list:
        df = df[df.crc.isin(crc)]

    if type(crc_rem) == list:
        df = df[~df.crc.isin(crc_rem)]

    if type(wap) == list:
        df = df[df.wap.isin(wap)]

    if type(wap_rem) == list:
        df = df[~df.wap.isin(wap_rem)]

    if type(take_type) == list:
        df = df[df.take_type.isin(take_type)]

    if type(use_type) == list:
        df = df[df.use_type.isin(use_type)]

    if type(catch_num) == list:
        df = df[df.catch_grp.isin(catch_num)]

    if type(catch_name) == list:
        df = df[df.catch_name.isin(catch_name)]

    if type(catch_grp_name) == list:
        df = df[df.catch_grp_.isin(catch_grp_name)]

    if type(swaz) == list:
        df = df[df.swaz.isin(swaz)]

    if type(gwaz) == list:
        df = df[df.gwaz.isin(gwaz)]

    if type(swaz_grp) == list:
        df = df[df.swaz_grp.isin(swaz_grp)]

    if type(cwms_zone) == list:
        df = df[df.cwms.isin(cwms_zone)]

    if type(gr_than) == list:
        df = df[df[allo_col[0]] > (gr_than[0] * 60*60*24*212/1000)]

    if type(years) == list:
        yrs_index = df.date.astype('str').str[0:4].astype('int')
        df = df[yrs_index.isin(years)]

    if isinstance(from_date, str):
        df = df[df.date >= from_date]

    if isinstance(to_date, str):
        df = df[df.date <= to_date]

    if sd_only:
        df = df[((df.sd1_150 > 0) & (df.take_type == 'Take Groundwater')) | (df.take_type == 'Take Surface Water')]

    df_grp = df.groupby(grp_by)
    df1 = df[df[use_col].notnull().values]
    df1_grp = df1.groupby(grp_by)

    df_sum = df_grp[allo_col].sum()
    df_count = df_grp.crc.nunique()
    df1_count = df1_grp.crc.nunique()

    df_col = ['tot_' + x for x in allo_col]
    df_sum.columns = df_col
    df1_sum = df1_grp[use_col + allo_col].sum()

    for i in range(len(allo_col)):
        df1_sum['usage/' + allo_col[i]] = df1_sum[use_col].div(df1_sum[allo_col[i]], axis=0).round(2)

    df2 = concat([df1_sum, df_sum], axis=1)
#    df12.columns = ['usage_m3', 'usage_allo_m3', 'usage_ratio', 'tot_allo_m3']
    for i in range(len(allo_col)):
        df2['allo/' + df_col[0]] = df2[allo_col[i]].div(df2[df_col[i]], axis=0).round(2)

    df2['crc_with_use_count'] = df1_count
    df2['tot_crc_count'] = df_count

    if debug:
        if export:
            df.to_csv(export_path, index=False)
        return(df)
    else:
        if export:
            df2.to_csv(export_path)
        return(df2)
