# -*- coding: utf-8 -*-
"""
Functions for allocation and usage processing.
"""

##################################################
#### Data processing functions

status_codes = ['Terminated - Replaced', 'Issued - Active', 'Terminated - Surrendered', 'Terminated - Cancelled', 'Terminated - Expired', 'Terminated - Lapsed', 'Issued - s124 Continuance']

allo_gis_dict = {'crc': 'The consent number.', 'take_type': 'The take type of the consent (i.e. Take Groundwater or Take Surface Water).', 'allo_block': 'The allocation block.', 'wap': 'The Water Abstraction Point (either Well or SWAP).', 'from_month': 'The month of the year that the consent can begin to take.', 'in_gw_allo': 'Should the consent be included in the GW allocation? Essentially, if the GW take is consumptive and not duplicated.', 'in_sw_allo': 'A combination of if a GW take is stream depleting and should be included in the SW allocation, or if a SW take is consumptive and not duplicated.', 'irr_area': 'The irrigation area (if known) for the irrigation take in hectares).', 'max_rate': 'The pro-rata max rate for the WAPs of the consent. The max rate of the overall consent distributed over the WAPs weighted by the individual WAPs max rate (in l/s).', 'max_rate_wap': 'The individual WAPs max rates (in l/s)', 'max_vol': 'The pro-rata max volume (in m3). Similar distribution as max_rate.', 'return_period': 'The number of days that the max_volume covers.', 'sd1_150': 'The stream depletion percent of the max rate.', 'use_type': 'The land use type that the consent/WAP is used for.', 'from_date': 'The start date of the consent.', 'to_date': 'The end date of the consent.', 'status_details': 'The specific and current status of the consent.', 'max_rate_crc': 'The max rate of the consent.', 'max_vol_crc': 'The max volume of the consent (in m3).', 'cav_crc': 'The consented annual volume (in m3).', 'min_flow': 'Does the consent have a min flow requirement?', 'cav': 'The estimated consented annual volume if a consent does not already have one (in m3). First, the estimate will be from the max_volume if one exists, if not then it is estimated from the max_rate.', 'daily_vol': 'An estimated daily volume. First from the max_volue if it exists, otherwise from the max_rate (in m3).'}


### Allocation and usage


def allo_proc(in_allo=True, corr_csv='S:/Surface Water/shared/base_data/database/database_corrections.csv', use_type_csv='S:/Surface Water/shared/base_data/usage/use_type_def.csv', export_path=None):
    """
    Function to take the consent allocation data from various sources, clean up the data, and reorganize the results to be used in the time series processing script. The most variety that can occur with one consent is multiple take types with multiple use types.

    Arguments:\n
    corr_csv -- csv with data corrections.\n
    use_type_csv -- csv where the use types are defined.\n
    export_path -- The path and csv file name that should be exported.
    """

    from numpy import in1d, nan, ceil
    from pandas import merge, to_datetime, read_csv, to_numeric, concat
#    from core.allo_use import wqn10
    from core.ecan_io import rd_sql
    from core.misc import save_df
    from core.ecan_io.SQL_databases import sql_arg

    #### Parameters
    take_types = ['Take Groundwater', 'Take Surface Water']
    sql1 = sql_arg()

    #### Load RAW data
#    takes = rd_sql(code='crc_acc')
    crc_wap = rd_sql(where_col='Activity', where_val=take_types, **sql1.get_dict('crc_wap_allo'))
    crc_cav = rd_sql(where_col='Activity', where_val=take_types, **sql1.get_dict('crc_cav'))
    crc_use_type = rd_sql(where_col='Activity', where_val=take_types, **sql1.get_dict('crc_use_type'))
    sd = rd_sql(**sql1.get_dict('sd'))

    dates = rd_sql(**sql1.get_dict('crc_dates'))[['crc', 'from_date', 'to_date', 'status_details']]

#    irr_dem = rd_sql(code='irr_dem_par')

    crc_rel = rd_sql(**sql1.get_dict('crc_relation'))

    use_type_code = read_csv(use_type_csv)
    val3 = [use_type_code.loc[i, ['values']].values[0].split(', ') for i in use_type_code.index]
    use_type_dict = dict(zip(use_type_code.use_type_code.tolist(), val3))

    #### Correct data errors
#    dict1 = {'crc_wap_act_acc': crc_wap, 'crc_details_acc': dates, 'crc_relation': crc_rel, 'crc_use_type_acc': crc_use_type, 'crc_gen_acc': crc_cav}
#
#    data_corr = read_csv(corr_csv)
#
#    for i in data_corr.index:
#        t1 = dict1[data_corr.Code[i]]
#        t1.loc[t1[data_corr.where_field[i]] == data_corr.where_field_val[i], data_corr.corr_field[i]] = to_datetime(data_corr.corr_field_val[i], errors='ignore')
#        dict1[data_corr.Code[i]] = t1
#
#    #### Reassign objects
#    crc_wap = dict1['crc_wap_act_acc'].copy()
#    crc_cav = dict1['crc_gen_acc'].copy()
#    crc_use_type = dict1['crc_use_type_acc'].copy()
#    dates = dict1['crc_details_acc'].copy()
##    irr_dem = dict1['irr_dem_par']
#    crc_rel = dict1['crc_relation'].copy()

    #### Initial processing post data correction

    ### crc_wap
    crc_wap.loc[:, 'in_sw_allo'] = crc_wap.loc[:, 'in_sw_allo'].str.upper()
#    crc_wap = crc_wap[(crc_wap.take_type == 'Take Groundwater') | ((crc_wap.take_type == 'Take Surface Water') & (crc_wap.in_sw_allo == 'YES'))].drop('in_sw_allo', axis=1)
    crc_wap.loc[crc_wap.wap == 'Migration: Not Classified', 'wap'] = nan
    crc_wap.loc[:, 'wap'] = crc_wap.loc[:, 'wap'].str.upper()
    crc_wap.loc[:, 'wap'] = crc_wap.wap.str.strip()
    crc_wap.loc[crc_wap.allo_block == 'Migration: Not Classified', 'allo_block'] = nan
    crc_wap.loc[:, 'max_rate'] = to_numeric(crc_wap.loc[:, 'max_rate'], errors='coerce')
    crc_wap.loc[:, 'max_rate_wap'] = to_numeric(crc_wap.loc[:, 'max_rate_wap'], errors='coerce')
    crc_wap.loc[:, 'max_vol'] = to_numeric(crc_wap.loc[:, 'max_vol'], errors='coerce')
    crc_wap.loc[:, 'sd'] = to_numeric(crc_wap.loc[:, 'sd'], errors='coerce')
    crc_wap.loc[:, 'return_period'] = to_numeric(crc_wap.loc[:, 'return_period'], errors='coerce')

    crc_wap.loc[:, 'from_month'] = crc_wap.loc[:, 'from_month'].str.upper()
    mon_index1 = (crc_wap.from_month == 'JUL') | (crc_wap.from_month == 'OCT')
    crc_wap.loc[~mon_index1, 'from_month'] = 'JUL'

    crc_wap.loc[crc_wap.max_vol <= 0, 'max_vol'] = nan
    crc_wap.loc[crc_wap.max_rate_wap <= 0, 'max_rate_wap'] = nan
    crc_wap.loc[crc_wap.max_rate <= 0, 'max_rate'] = nan
    crc_wap.loc[crc_wap.return_period <= 0, 'return_period'] = nan
    crc_wap.loc[crc_wap.max_rate.isnull(), 'max_rate'] = crc_wap.loc[crc_wap.max_rate.isnull(), 'max_rate_wap']

    ### crc_use_type
    crc_use_type.loc[:, 'in_gw_allo'] = crc_use_type.loc[:, 'in_gw_allo'].str.upper()
    crc_use_type.loc[crc_use_type.allo_block == 'Migration: Not Classified', 'allo_block'] = nan

    ## Redefine use types
    for i in use_type_dict:
        index1 = in1d(crc_use_type.use_type, use_type_dict[i])
        crc_use_type.loc[index1, 'use_type'] = i
    crc_use_type.loc[~in1d(crc_use_type.use_type, use_type_dict.keys()), 'use_type'] = 'other'

    ### Make sure allocation blocks are consistent between the two tables
    crc_wap_grp1 = crc_wap.groupby(['crc', 'take_type'])['allo_block']
    use_type_grp1 = crc_use_type.groupby(['crc', 'take_type'])['allo_block']

    count1 = crc_wap_grp1.count()
    count1.name = 'wap'
    count2 = use_type_grp1.count()
    count2.name = 'use'
    counts = concat([count1, count2], axis=1).dropna()
    counts_index1 = counts[(counts.use == 1)].index.get_values()
    counts_index2 = counts[(counts.use > 1) & (counts.use == counts.wap)].index.get_values()

    use_type0 = crc_use_type.set_index(['crc', 'take_type'])
    crc_wap0 = crc_wap.set_index(['crc', 'take_type'])

    crc_wap1 = crc_wap0.loc[counts_index1]
    crc_wap0.loc[counts_index1, 'allo_block'] = use_type0.loc[crc_wap1.index.get_values(), 'allo_block']
    crc_wap0.loc[counts_index2, 'allo_block'] = use_type0.loc[counts_index2, 'allo_block']

    ### Combine the use types with the WAPs
    use_type0 = use_type0.set_index(['allo_block'], append=True)
    crc_wap0 = crc_wap0.set_index(['allo_block'], append=True)
    crc_wap_grp = crc_wap.groupby(['crc', 'take_type', 'allo_block'])
    n_wap_set = crc_wap_grp['row_index'].count()
    n_wap_set.name = 'n_wap_set'
    use_type_grp = crc_use_type.groupby(['crc', 'take_type', 'allo_block'])
    n_set = use_type_grp['row_index'].count()
    n_set.name = 'n_set'

    ## One use type per crc
    one1 = n_set[n_set == 1].index.get_values()
    one2 = use_type0.loc[one1, ['in_gw_allo', 'use_type', 'irr_area']]
    crc_wap1 = crc_wap0.loc[one1]
    crc_wap1.index.set_names(crc_wap0.index.names, inplace=True)
    index2 = crc_wap1.index.get_values()

    one3 = concat([crc_wap1, one2.loc[index2]], axis=1)

    ## Same multi to multi
    both1 = concat([n_wap_set, n_set], axis=1)
    both2 = both1[(both1.n_set > 1) & (both1.n_set == both1.n_wap_set)].index.get_values()

    both3 = concat([crc_wap0.loc[both2], use_type0.loc[both2, ['in_gw_allo', 'use_type', 'irr_area']]], axis=1)

    ## Unequal sets
    dis1 = both1[(both1.n_set > 1) & (both1.n_set != both1.n_wap_set)]
    dis2 = dis1[dis1.n_wap_set.notnull()]

    # Agg WAPs
    dis3 = crc_wap0.loc[dis2.index.get_values()].set_index('wap', append=True)
    dis_grp = dis3.reset_index().sort_values(['crc', 'row_index']).groupby(['crc', 'take_type', 'allo_block', 'wap'])
    sum1 = dis_grp['max_rate_wap', 'max_rate', 'max_vol'].sum()
    mean1 = dis_grp['return_period', 'sd'].mean()
    first1 = dis_grp['in_sw_allo', 'from_month'].first()
    new1 = concat([sum1, mean1, first1], axis=1).reset_index('wap')
    new_index1 = new1.groupby(level=['crc', 'take_type', 'allo_block']).first().index.get_values()
    new_use1 = use_type0.loc[new_index1]
    new_use_count = new_use1.reset_index().sort_values(['crc', 'row_index']).groupby(['crc', 'take_type', 'allo_block'])['row_index'].count()
    new_use_count.name = 'new_use'
    new_wap_count = new1.groupby(level=['crc', 'take_type', 'allo_block'])['from_month'].count()
    new_wap_count.name = 'new_wap'
    new_both1 = concat([new_wap_count, new_use_count], axis=1)
    new_both2 = new_both1[new_both1.new_wap == new_both1.new_use].index.get_values()

    new_both3 = concat([new1.loc[new_both2], new_use1.loc[new_both2]], axis=1)

    # Agg use_types - and remove use_types that are not 'in_gw_allo'
    mis1 = new_both1[new_both1.new_wap != new_both1.new_use].index.get_values()
    last_use1 = use_type0.loc[mis1].drop('row_index', axis=1)
    last_use2 = last_use1.set_index('use_type', append=True)
    last_use3 = last_use2[(last_use2.index.get_level_values('take_type') == 'Take Surface Water') | (last_use2.in_gw_allo == 'YES')]
    last_use4 = last_use3.groupby(level=['crc', 'take_type', 'allo_block', 'use_type'])['irr_area'].sum()
    last_use4 = concat([last_use4, last_use3.groupby(level=['crc', 'take_type', 'allo_block', 'use_type'])['in_gw_allo'].first()], axis=1)
    last_use4 = last_use4.reset_index('use_type')
    last_count1 = last_use4.groupby(level=['crc', 'take_type', 'allo_block'])['use_type'].count()
    last_count2 = last_count1[last_count1 == 1].index.get_values()
    new3 = new1.loc[last_count2]

    new_agg_use = concat([new3, last_use4.loc[new3.index.get_values()]], axis=1)

    # Randomly distribute the rest
    last_count3 = last_count1[last_count1 > 1]
    last_use5 = last_use4.loc[last_count3.index.get_values()]
    last_wap1 = new1.loc[last_count3.index.get_values()]

    new_lst = []
    for i in last_count3.index.get_values():
        last_wap1 = new1.loc[i]
        last_use5 = last_use4.loc[i]
        len_wap = len(last_wap1)
        len_use = len(last_use5)
        ceil1  = int(ceil(len_wap/float(len_use)))
        diff1 = len_use * ceil1 - len_wap

        new_use = concat([last_use5] * ceil1)
        new_use = new_use[:(len(new_use) - diff1)]

        new_both1 = concat([last_wap1, new_use], axis=1)
        new_lst.append(new_both1)

    new_both2 = concat(new_lst)

    ## Combine all data together
    crc_wap_use = concat([one3.drop('row_index', axis=1), both3.drop('row_index', axis=1), new_both3.drop('row_index', axis=1), new_agg_use, new_both2]).reset_index()

    ## Find the missing WAPs per consent
    crc_wap_mis1 = crc_wap_use.loc[crc_wap_use.wap.isnull(), 'crc'].unique()
    crc_wap4 = crc_wap_use[['crc', 'wap']]

    for i in crc_wap_mis1:
        crc1 = crc_rel[in1d(crc_rel.crc, i)].crc_child.values
        wap1 = []
        while (len(crc1) > 0) & (len(wap1) == 0):
            wap1 = crc_wap4.loc[in1d(crc_wap4.crc, crc1), 'wap'].values
            crc1 = crc_rel[in1d(crc_rel.crc, crc1)].crc_child.values
        if len(wap1) > 0:
            crc_wap_use.loc[crc_wap_use.crc == i, 'wap'] = wap1[0]

    ### dates - only use crc's with from and to dates
    dates.loc[:, 'from_date'] = to_datetime(dates.from_date, errors='coerce')
    dates.loc[:, 'to_date'] = to_datetime(dates.to_date, errors='coerce')
#    dates = dates.dropna()

    ## Append the dates
    crc_wap_use_dates = merge(crc_wap_use, dates, on='crc', how='left')

    ### CAV processing
    crc_cav.loc[:, 'cav'] = to_numeric(crc_cav.loc[:, 'cav'], errors='coerce')
    crc_cav.loc[:, 'min_flow'] = crc_cav.loc[:, 'min_flow'].str.upper()
    crc_cav.loc[crc_cav.cav < 200, 'cav'] = nan
    crc_cav.columns = ['crc', 'take_type', 'cav_crc', 'min_flow']

    sum_wap = crc_wap_use_dates.groupby(['crc', 'take_type'])['max_rate_wap', 'max_rate', 'max_vol'].sum().round(1)
    sum_wap.columns = ['max_rate_wap_crc', 'max_rate_crc', 'max_vol_crc']
    cav0 = crc_wap_use_dates.set_index(['crc', 'take_type']).loc[sum_wap.index.get_values()]
    cav1 = concat([cav0, sum_wap.loc[cav0.index.get_values()]], axis=1).reset_index()
    cav1['rate_ratio'] = cav1['max_rate'] / cav1['max_rate_crc']

    ## Combine cav with other data
    cav2 = merge(cav1, crc_cav, on=['crc', 'take_type'], how='left')
    cav2['cav'] = (cav2['cav_crc'] * cav2['rate_ratio']).round(2)

    ## Estimate other cav's and daily vol
    cav2['daily_vol'] = (cav2['max_vol'] / cav2['return_period']).round(2)
    cav2.loc[cav2.daily_vol.isnull(), 'daily_vol'] = (cav2.loc[cav2.daily_vol.isnull(), 'max_rate'] * 60*60*24/1000).round(2)
    oct_index = cav2.cav.isnull() & (cav2.from_month == 'OCT')
    cav2.loc[oct_index, 'cav'] = cav2.loc[oct_index, 'daily_vol'] * 212
    cav2.loc[cav2.cav.isnull(), 'cav'] = cav2.loc[cav2.cav.isnull(), 'daily_vol'] * 365
    allo = cav2.drop('rate_ratio', axis=1).drop_duplicates(subset=['crc', 'take_type', 'allo_block', 'wap', 'use_type'])

    #### Add in the stream depletion percentages
    allo2 = merge(allo.drop('sd', axis=1), sd[['wap', 'sd1_150']], on='wap', how='left')

    #### Reorder the columns
    cols = allo2.columns.tolist()
    init_cols = ['crc', 'take_type', 'allo_block', 'wap', 'use_type', 'max_rate', 'daily_vol', 'cav', 'max_vol', 'return_period', 'from_date', 'to_date', 'status_details']
    t1 = [cols.remove(i) for i in init_cols]
    cols[0:0] = init_cols
    allo3 = allo2[cols].copy()

    #### Calc irr area where not populated
    irr_crc1 = allo3[allo3.use_type == 'irrigation'].copy()
    irr_crc2 = irr_crc1.groupby(['crc'])[['irr_area', 'daily_vol']].sum()
    irr_crc3 = irr_crc2[irr_crc2.irr_area > 0].copy()
    irr_crc3.loc[:, 'irr_day_ratio'] = irr_crc3['daily_vol']/(irr_crc3.irr_area * 10)
    irr_crc4 = irr_crc3[(irr_crc3['irr_day_ratio'] >= 2) & (irr_crc3['irr_day_ratio'] <= 6)]
    mm_hect = irr_crc4['irr_day_ratio'].mean()

    allo3.loc[:, 'irr_area_est'] = allo3['irr_area']
    allo3.loc[~allo3.crc.isin(irr_crc4.index), 'irr_area_est'] = (allo3.loc[~allo3.crc.isin(irr_crc4.index), 'daily_vol'] / mm_hect / 10).round(1)

    #### Make corrections
    allo3.loc[(allo3['take_type'] == 'Take Groundwater') & (allo3['in_sw_allo'] == 'YES'), 'in_gw_allo'] = 'YES'

    #### Save data and return object
    if isinstance(export_path, str):
        save_df(allo3, export_path, index=False)

    return(allo3)


def allo_gis_proc(allo, export=True, export_shp='allo_gis.shp', export_csv='allo_gis.csv'):
    """
    Function to assign locations to the WAPs and to determine the GIS areas of each WAP.
    """
    from core.ecan_io import rd_sql
    from core.spatial import xy_to_gpd
    from numpy import in1d, arange
    from geopandas.tools import sjoin
    from pandas import merge, concat
    from core.ecan_io.SQL_databases import sql_arg

    #### Read data
    sql1 = sql_arg()
    swaz = rd_sql(**sql1.get_dict('swaz_gis'))
    gwaz = rd_sql(**sql1.get_dict('gwaz_gis'))
    catch = rd_sql(**sql1.get_dict('catch_gis'))
    cwms = rd_sql(**sql1.get_dict('cwms_gis'))
    poly_dict = {'swaz': swaz, 'gwaz': gwaz, 'catch': catch, 'cwms': cwms}

    wap_loc1 = rd_sql(**sql1.get_dict('well_details'))[['wap', 'NZTMX', 'NZTMY']]
    wap_loc = xy_to_gpd('wap', 'NZTMX', 'NZTMY', wap_loc1)
    crc_loc = rd_sql(**sql1.get_dict('crc_gis'))[['crc', 'geometry']]

    #### Prepare allo data
#    allo1 = allo[in1d(allo.status_details, status_codes)]
    mis_waps = allo[allo.wap.isnull()]
    all_waps = allo.wap.unique()

    #### Assign locations
    all_wap_loc = wap_loc[in1d(wap_loc.wap, all_waps)]
    mis_wap_loc = crc_loc[in1d(crc_loc.crc, mis_waps.crc.unique())]

    #### Create new WAP numbers
    x1 = mis_wap_loc.geometry.apply(lambda x: x.x).round(3)
    x1.name = 'x'
    y1 = mis_wap_loc.geometry.apply(lambda x: x.y).round(3)
    y1.name = 'y'
    mis_wap_xy = concat([mis_wap_loc['crc'], x1, y1], axis=1)
    t1 = mis_wap_xy[['x', 'y']].drop_duplicates()
    t1.loc[:, 'wap_new'] = arange(len(t1)) + 10000
    mis_wap_xy1 = merge(mis_wap_xy, t1, on=['x', 'y'], how='left').drop(['x', 'y'], axis=1)
    mis_wap_loc2 = merge(mis_wap_loc, mis_wap_xy1, on='crc', how='left')

    #### add areas to locations
    all_wap_loc1 = all_wap_loc.copy()
    mis_wap_loc1 = mis_wap_loc2.copy()
    for i in poly_dict:
        poly = poly_dict[i]
        all_wap_loc1 = sjoin(all_wap_loc1, poly, how='left', op='within').drop('index_right', axis=1)
        mis_wap_loc1 = sjoin(mis_wap_loc1, poly, how='left', op='within').drop('index_right', axis=1)

    #### Put GIS data into the allo object
    allo1 = allo.copy()
#    allo1 = merge(allo1, all_wap_loc1.drop('geometry', axis=1), on='wap', how='left')
#    allo2 = merge(allo, mis_wap_loc1.drop('geometry', axis=1), on='crc')
#    allo2.loc[allo2.wap.isnull(), 'wap'] = allo2.loc[allo2.wap.isnull(), 'wap_new']
#    allo3 = concat([allo1[~in1d(allo1.crc, mis_wap_loc1.crc)], allo2]).drop('wap_new', axis=1)

    gis1 = merge(all_wap_loc1, allo1, on='wap', how='left')
    gis2 = merge(mis_wap_loc1, allo1, on='crc', how='left')
    gis2.loc[gis2.wap.isnull(), 'wap'] = gis2.loc[gis2.wap.isnull(), 'wap_new']

    gis3 = concat([gis1, gis2]).drop('wap_new', axis=1)

    gis3.loc[:, ['from_date', 'to_date']] = gis3.loc[:, ['from_date', 'to_date']].astype('str')

    ## Add in the x and y coordinates
    gis3['x'] = gis3.geometry.apply(lambda i: i.x)
    gis3['y'] = gis3.geometry.apply(lambda i: i.y)

    ## Reorder columns
    cols = gis3.columns.tolist()
    init_cols = ['crc', 'take_type', 'allo_block', 'wap', 'use_type', 'max_rate', 'daily_vol', 'cav', 'max_vol', 'return_period', 'from_date', 'to_date', 'status_details']
    t1 = [cols.remove(i) for i in init_cols]
    cols[0:0] = init_cols
    gis4 = gis3[cols]

    #### Save data
    if export:
        gis4.to_file(export_shp)
        gis4.drop('geometry', axis=1).to_csv(export_csv, encoding='utf-8', index=False)
    return(gis4)


def allo_ts_apply(wap, start_date='2014-07-01', end_date='2016-06-30', from_col='from_date', to_col='to_date', freq='D', mon_col='from_month', daily_vol_col='daily_vol', cav_col='cav'):
    """
    Function that converts the allocation data to a monthly time series.
    """
    from pandas import DataFrame, DateOffset, date_range, Timestamp, Series
    from numpy import nan, in1d

    from_date = Timestamp(wap[from_col])
    to_date = Timestamp(wap[to_col])
    mon = wap[mon_col]
    start = Timestamp(start_date)
    end = Timestamp(end_date)
    take_type = wap.name[1]

    if from_date < start:
        from_date = start
    if to_date > end:
        to_date = end

    if (freq == 'A') & (mon == 'OCT'):
        dates1 = date_range(from_date, to_date - DateOffset(2) + DateOffset(months=9), freq='A-APR')
    elif (freq == 'A') & (mon == 'JUL'):
        dates1 = date_range(from_date, to_date - DateOffset(2) + DateOffset(years=1), freq='A-JUN')
    elif freq == 'sw_rates':
        dates1 = date_range(from_date, to_date, freq='AS-JAN')
    else:
        dates1 = date_range(from_date, to_date, freq=freq)

    if freq == 'sw_rates':
        if len(dates1) > 0:
            if (take_type == 'Take Surface Water'):
                s1 = Series(wap['max_rate'], index=dates1)
                return(s1)
            elif ((take_type == 'Take Groundwater') & (wap['min_flow'] == 'YES')):
                s1 = Series(wap['max_rate'] * wap['sd1_150'] * 0.01, index=dates1).round(2)
                return(s1)

    elif mon == 'OCT':
        dates_index = in1d(dates1.month, [10, 11, 12, 1, 2, 3, 4])
        if len(dates1) > 0:
            if freq == 'D':
                s1 = Series(0, index=dates1)
                s1.loc[dates_index] = wap[daily_vol_col]
            elif freq == 'M':
                start_n_days = (dates1[0] - from_date).days + 1
                if start_n_days > dates1[0].day:
                    start_n_days = dates1[0].day
                end_n_days = dates1[-1].day - (dates1[-1] - to_date).days
                if end_n_days > dates1[-1].day:
                    end_n_days = dates1[-1].day
                days = dates1.days_in_month.values
                days[0] = start_n_days
                days[-1] = end_n_days
                vol = (days/213.0 * wap[cav_col]).round(1)
                s1 = Series(0, index=dates1)
                s1.loc[dates_index] = vol[dates_index]
            elif freq == 'A':
                extra_days = 91.0
                start_n_days = (dates1[0] - from_date).days
                if start_n_days > (dates1[0].dayofyear + extra_days):
                    start_n_days = (dates1[0].dayofyear + extra_days)
                end_n_days = (dates1[-1].dayofyear + extra_days) - (dates1[-1] - to_date).days
                if end_n_days > (dates1[-1].dayofyear + extra_days):
                    end_n_days = (dates1[-1].dayofyear + extra_days)
                dayofyear = dates1.dayofyear.values
                days = dayofyear + extra_days
                if len(days) == 1:
                    days[0] = days[0] - (days[0] - start_n_days) - (days[0] - end_n_days)
                else:
                    days[0] = start_n_days
                    days[-1] = end_n_days
                vol = (days/(dayofyear + extra_days) * wap[cav_col]).round(1)
                s1 = Series(vol, index=dates1 + DateOffset(months=2))
            return(s1)
    elif mon == 'JUL':
        if freq == 'D':
            s1 = Series(wap[daily_vol_col], index=dates1)
        elif freq == 'M':
            start_n_days = (dates1[0] - from_date).days
            end_n_days = dates1[-1].day - (dates1[-1] - to_date).days
            days = dates1.days_in_month.values
            days[0] = start_n_days
            days[-1] = end_n_days
            vol = (days/365.0 * wap[cav_col]).round(1)
            s1 = Series(vol, index=dates1)
        elif freq == 'A':
            extra_days = 183.0
            start_n_days = (dates1[0] - from_date).days
            end_n_days = (dates1[-1].dayofyear + extra_days) - (dates1[-1] - to_date).days
            dayofyear = dates1.dayofyear.values
            days = dayofyear + extra_days
            if len(days) == 1:
                days[0] = days[0] - (days[0] - start_n_days) - (days[0] - end_n_days)
            else:
                days[0] = start_n_days
                days[-1] = end_n_days
            vol = (days/(dayofyear + extra_days) * wap[cav_col]).round(1)
            s1 = Series(vol, index=dates1)
        return(s1)


def allo_filter(allo, start='1900-07-01', end='2020-06-30', from_col='from_date', to_col='to_date', in_allo=True):
    """
    Function to take an allo DataFrame and filter out the consents that cannot be converted to a time series due to missing data.
    """
    from pandas import to_datetime, Timestamp

    allo1 = allo.copy()
    allo1.loc[:, to_col] = to_datetime(allo1.loc[:, to_col])
    allo1.loc[:, from_col] = to_datetime(allo1.loc[:, from_col])

    ### Remove consents without daily volumes (and consequently yearly volumes)
    allo2 = allo1[allo1.daily_vol.notnull()]

    ### Remove consents without to/from dates or date ranges of less than a month
    allo3 = allo2[allo2.from_date.notnull() & allo2.to_date.notnull()]

    ### Restrict dates
    start_time = Timestamp(start)
    end_time = Timestamp(end)

    allo4 = allo3[(allo3.to_date - start_time).dt.days > 31]
    allo5 = allo4[(end_time - allo4.from_date).dt.days > 31]

    allo5 = allo5[(allo5.to_date - allo5.from_date).dt.days > 31]

    ### Restrict by status_details
    allo6 = allo5[allo5.status_details.isin(['Terminated - Replaced', 'Issued - Active', 'Terminated - Surrendered', 'Terminated - Expired', 'Terminated - Lapsed', 'Issued - s124 Continuance', 'Terminated - Cancelled'])]

    ### Remove Hydroelectric consents
    allo6 = allo6[allo6.use_type != 'hydroelectric']

    ### In allocation columns
    if in_allo:
        allo6 = allo6[(allo6.take_type == 'Take Surface Water') | ((allo6.take_type == 'Take Groundwater') & (allo6.in_gw_allo == 'YES'))]
        allo6 = allo6[(allo6.take_type == 'Take Groundwater') | ((allo6.take_type == 'Take Surface Water') & (allo6.in_sw_allo == 'YES'))]

    ### Index the DataFrame
    allo7 = allo6.set_index(['crc', 'take_type', 'allo_block', 'wap'])

    ### Return
    return(allo7)


def allo_ts_proc(allo, start='1900-07-01', end='2020-06-30', in_allo=True, from_col='from_date', to_col='to_date', freq='M', mon_col='from_month', daily_vol_col='daily_vol', cav_col='cav', export_path=None):
    """
    Combo function to completely create a time series from the allocation DataFrame.
    """
    from core.allo_use import allo_ts_apply, allo_filter
    from core.misc import save_df

    allo1 = allo.copy()
    allo1 = allo1.drop_duplicates(subset=['crc', 'take_type', 'allo_block', 'wap'])

    allo2 = allo_filter(allo1, start=start, end=end, in_allo=in_allo)

    allo3 = allo2.apply(lambda x: allo_ts_apply(x, start_date=start, end_date=end, freq=freq, from_col=from_col, to_col=to_col, mon_col=mon_col, daily_vol_col=daily_vol_col, cav_col=cav_col), axis=1)

    allo4 = allo3.stack().reset_index()
    allo4.columns = ['crc', 'take_type', 'allo_block', 'wap', 'date', 'allo']

    if isinstance(export_path, str):
        save_df(allo4, export_path, index=False)

    return(allo4)


def wqn10(irr1):
    """
    Official Ecan method to estimate yearly water demand for water allocation.
    """
    irr = irr1.copy().drop_duplicates()
    irr.loc[:, 'ann1'] = 910 - 1.6*(irr.paw - 100)
    irr.loc[irr.paw < 100, 'ann1'] = 750
    irr.loc[irr.paw > 200, 'ann1'] = 910

    irr.loc[:, 'demand'] = irr.ann1 - irr.siwl
    irr.loc[:, 'ann_vol_wqn'] = irr.demand * 10 * irr.irr_area
    irr.loc[irr.ann_vol_wqn.isnull(), 'ann_vol_wqn'] = 0

    irr2 = irr[['crc', 'ann_vol_wqn']].groupby('crc').max().reset_index()
    irr2 = irr2[irr2.ann_vol_wqn > 100]

    return(irr2)

### Usage


def allo_use_proc(allo_ts_mon, usage, export_path=None):
    """
    Function to merge the usage data with the allocation data.
    """
    from pandas import merge, to_datetime
    from core.ts import grp_ts_agg
    from core.misc import save_df

    #### Process the usage data
    use1 = grp_ts_agg(usage, 'wap', 'date', 'M').sum().reset_index()
#    crc_use1 = use1[use1.wap.str.contains('CRC')]
    use1.loc[:, 'usage'] = use1.loc[:, 'usage'].round(2)
#    use1.loc[:, 'date'] = use1.loc[:, 'date'].astype(str)

    #### Process allo data for the distribution of usage data
    allo_ts = allo_ts_mon.copy()
    allo_ts.loc[:, 'date'] = to_datetime(allo_ts.loc[:, 'date'])

    ### multiple consents on the same month with the same WAP
    allo_ts['dup_wap']  = allo_ts.groupby(['wap', 'date'])['allo'].transform('count')

    ### Combine use with allo/wap
    allo_use1 = merge(allo_ts, use1, on=['wap', 'date'], how='left')
    use2 = (allo_use1.loc[:, 'usage'] / allo_use1.loc[:, 'dup_wap']).round(2)
    allo_use1.loc[:, 'usage'] = use2

    ### Remove unnecessary column
    allo_use2 = allo_use1.drop('dup_wap', axis=1)

    #### Export and return
    if isinstance(export_path, str):
        save_df(allo_use2, export_path, index=False)

    return(allo_use2)


def est_use(allo_use, allo_use_ros, allo_gis, date_col='date', usage_col='mon_usage_m3', allo_col='ann_restr_allo_m3', min_ratio=0.2, date_start='2012-06-01', date_end='2015-06-30', export=False, export_path='usage_allo_est1.h5'):
    """
    Function to estimate the usage during a period that has a decent amount of actual usage data.
    """
#    def resample1(df):
#        df.index = df.dates
#        df_grp = df.resample('A-JUN')
#        df2 = df_grp['mon_vol'].transform(sum).round(2)
#        return(df2)

    from pandas import merge, concat, date_range, to_datetime, TimeGrouper
    from numpy import nan
    from core.ts import grp_ts_agg

    ## Make sure dates are as datetime
    allo_use.loc[:, date_col] = to_datetime(allo_use[date_col])
    allo_use_ros.loc[:, date_col] = to_datetime(allo_use_ros[date_col])

    allo_use.rename(columns={'allo': 'mon_allo_m3', 'usage': 'mon_usage_m3'}, inplace=True)

    ## Resample data
    ann_allo1 = grp_ts_agg(allo_use, ['crc', 'take_type', 'allo_block', 'wap'], date_col, 'A-JUN').transform('sum')
    ann_allo1.rename(columns={'mon_allo_m3': 'ann_allo_m3', 'mon_usage_m3': 'ann_usage_m3'}, inplace=True)
    ann_allo = concat([allo_use[['crc', 'take_type', 'allo_block', 'wap', 'date']], ann_allo1.reset_index(drop=True)], axis=1)

    ## Merge usage and allo
    join_grp = ['crc', date_col, 'take_type', 'allo_block', 'wap']
    allo1 = merge(allo_use, ann_allo, on=join_grp, how='left')
    allo2 = merge(allo1, allo_use_ros.drop(['mon_allo_m3', 'mon_usage_m3'], axis=1), on=join_grp, how='left')

    allo3 = merge(allo2, allo_gis.drop(['max_rate', 'daily_vol', 'cav', 'cav_crc', 'from_date', 'to_date', 'max_rate_wap_crc', 'max_rate_crc', 'max_rate_wap', 'max_vol', 'max_vol_crc', 'return_period'], axis=1), on=['crc', 'take_type', 'allo_block', 'wap'], how='inner')
    allo4 = allo3.drop_duplicates(subset=join_grp)

    ## Combine annual allo to one column
    allo5 = allo4.copy()
    allo_index1 = allo5[allo_col].isnull()
    allo5.loc[allo_index1, allo_col] = allo5.loc[allo_index1, 'ann_allo_m3']
    allo5.loc[allo_index1, 'mon_restr_allo_m3'] = allo5.loc[allo_index1, 'mon_allo_m3']

    ## Estimate the usage ratios by catchment number
#    allo_use1['catch'] = floor(allo_use1.catchment_num * 0.0001)
#    allo_use1.loc[allo_use1.catch == 0, 'catch'] = nan
    allo5.rename(columns={'catch_grp':'catch'}, inplace=True)
    cols1 = ['crc', date_col, 'take_type', 'use_type', 'allo_block', 'wap', usage_col, allo_col, 'mon_restr_allo_m3', 'catch', 'cwms']
    allo6 = allo5[cols1].copy()
    allo6['usage_ratio_est'] = nan
    allo6['usage_est'] = nan

    ##  Index useful years
    dates = date_range(date_start, date_end, freq='M')
    dates_index = allo6[date_col].isin(dates)
    allo_use1a = allo6[dates_index]

    ## Estimate the usage ratio from annual updated allocations
    allo_use1b = allo_use1a[allo_use1a[usage_col].notnull()]

    ## Aggregate waps
    sum1 = allo_use1b.groupby(['crc', 'take_type', 'allo_block', date_col])[[usage_col, allo_col, 'mon_restr_allo_m3']].sum()
    first1 = allo_use1b.groupby(['crc', 'take_type', 'allo_block', date_col])[['use_type', 'catch', 'cwms']].first()
    allo_use2 = concat([sum1, first1], axis=1).reset_index()

    ## Calc usage ratio
    allo_use2.loc[:, 'usage_ratio'] = allo_use2[usage_col]/allo_use2[allo_col]
    allo_use2.loc[allo_use2.usage_ratio > 2 ,'usage_ratio'] = nan

    ## group by month, use type, and catchment number
    allo_use2_grp = allo_use2[[date_col, 'use_type', 'usage_ratio', 'catch', 'mon_restr_allo_m3']].groupby([date_col, 'use_type', 'catch'])

    allo_use3 = allo_use1a.copy()
    allo_use3_grp = allo_use3[[date_col, 'use_type', 'catch',  'mon_restr_allo_m3']].groupby([date_col, 'use_type', 'catch'])

    ## group by month, use type, and zone
    allo_use4_grp = allo_use2[[date_col, 'use_type', 'usage_ratio', 'cwms', 'mon_restr_allo_m3']].groupby([date_col, 'use_type', 'cwms'])
    allo_use5_grp = allo_use3[[date_col, 'use_type', 'cwms', 'mon_restr_allo_m3']].groupby([date_col, 'use_type', 'cwms'])

    ## group by month, use type, and region
    allo_use6_grp = allo_use2[[date_col, 'use_type', 'usage_ratio', 'mon_restr_allo_m3']].groupby([date_col, 'use_type'])
    allo_use7_grp = allo_use3[[date_col, 'use_type', 'mon_restr_allo_m3']].groupby([date_col, 'use_type'])

    ## Estimate mean usage ratio by groups

    # By catchment
    mean_use1 = allo_use2_grp['usage_ratio'].mean().round(4)
    count_use1 = allo_use2_grp['mon_restr_allo_m3'].count()

    count_tot1 = allo_use3_grp.count()

    count2 = concat([count_use1, count_tot1], axis=1)
    count2.columns = ['count1', 'count_tot']
    wap_ratio = count2.count1/count2.count_tot

    # By zone
    mean_use_zone = allo_use4_grp['usage_ratio'].mean().round(4)
    count_use_zone = allo_use4_grp['mon_restr_allo_m3'].count()

    count_tot_zone = allo_use5_grp.count()

    count_zone = concat([count_use_zone, count_tot_zone], axis=1)
    count_zone.columns = ['count1', 'count_tot']
    wap_ratio_zone = count_zone.count1/count_zone.count_tot

    # By region
    mean_use_cant = allo_use6_grp['usage_ratio'].mean().round(4)
    count_use_cant = allo_use6_grp['mon_restr_allo_m3'].count()

    count_tot_cant = allo_use7_grp.count()

    count_cant = concat([count_use_cant, count_tot_cant], axis=1)
    count_cant.columns = ['count1', 'count_tot']
    wap_ratio_cant = count_cant.count1/count_cant.count_tot

    ## Apply ratios and calc usage for all WAPs based on filter criteria
    usage_catch = concat([mean_use1, wap_ratio], axis=1)
    usage_catch.columns = ['usage_ratio_catch', 'wap_ratio_catch']
    usage_zone = concat([mean_use_zone, wap_ratio_zone], axis=1)
    usage_zone.columns = ['usage_ratio_zone', 'wap_ratio_zone']
    usage_cant = concat([mean_use_cant, wap_ratio_cant], axis=1)
    usage_cant.columns = ['usage_ratio_cant', 'wap_ratio_cant']

    allo1a = merge(allo_use1a, usage_catch.reset_index(), on=[date_col, 'use_type', 'catch'], how='left')
    allo2a = merge(allo1a, usage_zone.reset_index(), on=[date_col, 'use_type', 'cwms'], how='left')
    allo3 = merge(allo2a, usage_cant.reset_index(), on=[date_col, 'use_type'], how='left')

    no_usage_index = allo3[usage_col].isnull()
    min_catch_index = allo3['wap_ratio_catch'] >= min_ratio
    min_zone_index = allo3['wap_ratio_zone'] >= min_ratio

    allo3['usage_ratio_est'] = nan
    allo3.loc[(min_catch_index), 'usage_ratio_est'] = allo3.loc[(no_usage_index & min_catch_index), 'usage_ratio_catch']
    allo3.loc[(allo3.usage_ratio_est.isnull() & min_zone_index), 'usage_ratio_est'] = allo3.loc[(allo3.usage_ratio_est.isnull() & min_zone_index), 'usage_ratio_zone']
    allo3.loc[allo3.usage_ratio_est.isnull(), 'usage_ratio_est'] = allo3.loc[allo3.usage_ratio_est.isnull(), 'usage_ratio_cant']

    allo3['usage_est'] = allo3[usage_col].round(2)
    allo3.loc[no_usage_index, 'usage_est'] = (allo3.loc[no_usage_index, allo_col] * allo3.loc[no_usage_index, 'usage_ratio_est']).round(2)

    ### Merge all of the allo with the recent usage est
    allo5.loc[dates_index, 'usage_est'] = allo3['usage_est'].values
    allo5.loc[dates_index, 'usage_ratio_est'] = allo3['usage_ratio_est'].values

    col_names = allo2.columns.tolist()
#    col_names.extend(['sd1_150', 'usage_ratio_est', 'usage_est'])
    col_names.extend(['usage_ratio_est', 'usage_est'])

    allo_use2 = allo5[col_names].copy()

    if export:
        allo_use2.to_hdf(export_path, 'usage_est', mode='w')

    return(allo_use2)


def hist_sd_use(usage, allo_gis, vcn_grid2, et_col='pe', date_col='date', export=True, export_mon_path='sd_mon_vol.csv', export_reg_path='sd_reg.csv', export_sd_est_path='sd_est_all_mon_vol.csv'):
    """
    Function to create stream depletion for all consents from monthly sums of usage and ET.
    """
    from pandas import read_csv, to_datetime, datetime, merge, MultiIndex, DataFrame, date_range, IndexSlice, concat, Grouper, to_numeric
    from numpy import in1d, repeat, nan
    from core.stats import lin_reg
    from core.ecan_io.met import rd_niwa_vcsn
    from core.ts import grp_ts_agg
    from core.allo_use import allo_filter
    from core.misc import save_df

    ### Name parameters
    mon_allo_name = 'mon_allo_m3'
    ann_allo_name = 'ann_allo_m3'
    mon_restr_allo_name = 'mon_restr_allo_m3'
    ann_restr_allo_name = 'ann_restr_allo_m3'
    usage_name = 'usage_est'
    raw_usage_name = 'mon_usage_m3'
    sd_usage_name = 'sd_usage'
    unique_ids = ['crc', 'take_type', 'allo_block', 'wap']

    ### Read in data
    vcn_grid2.rename(columns={'niwa_id': 'site'}, inplace=True)
    vcn_sites = vcn_grid2[vcn_grid2.site.notnull()]
    vcn_sites_lst = vcn_sites.site.values
    vcn_data = rd_niwa_vcsn('PET', vcn_sites_lst, include_sites=True).drop(['x', 'y'], axis=1)
    vcn_data2 = merge(vcn_data, vcn_sites[['site', 'catch_grp']], on='site').drop('site', axis=1)

    usage_est = usage.copy()
    usage_est.loc[:, date_col] = to_datetime(usage_est[date_col])
    usage_est.rename(columns={date_col: 'time'}, inplace=True)


    ### Resample the ET series to monthly sums
    et_mon = grp_ts_agg(vcn_data2, 'catch_grp', 'time', 'M').mean().reset_index()
    et_mon.loc[:, 'catch_grp'] = to_numeric(et_mon.loc[:, 'catch_grp'], errors='coerce')

    ### Filter only the GW sites with SD and surface water sites
    allo_gis1 = allo_gis[allo_gis.sd1_150.notnull() | (allo_gis.take_type == 'Take Surface Water')].drop_duplicates(subset=unique_ids)
    allo_gis2 = allo_filter(allo_gis1).reset_index()
    usage1 = merge(usage_est, allo_gis2[unique_ids], on=unique_ids)

    first_date = et_mon.time[0]
    time_temp = usage1[usage1.usage_est.notnull()].time
    last_date = time_temp[time_temp.last_valid_index()]
    usage1 = usage1[(usage1.time <= last_date) & (usage1.time >= first_date)]

    ### Estimate the number of days per each month
    days1 = to_datetime(usage1.time.unique()).sort_values()
    days = days1[in1d(days1, et_mon.time.unique())]
#    count1 = days.days_in_month

    ### Merge allo/use with the spatial info
    tot_use1 = merge(usage1, allo_gis2, how='left')
    tot_use1 = tot_use1[tot_use1.catch_grp.notnull()]

    ### Calc SD usage
#    tot_use1['sd_usage'] = tot_use1['sd1_150'] * 0.01 * tot_use1['usage_est']
#    tot_use1.loc[(tot_use1.take_type == 'Take Surface Water'), 'sd_usage'] = tot_use1.loc[(tot_use1.take_type == 'Take Surface Water'), 'usage_est']

    ### Select ET data that coincide with the usage data
    et_mon_sel = et_mon[in1d(et_mon.catch_grp, tot_use1.catch_grp.unique()) & (et_mon.catch_grp != nan)]
    et_mon_sel_rec = et_mon_sel[et_mon_sel.time.isin(days)]

    ### Re-group the use types
    tot_use1['use_type_reg'] = tot_use1['use_type']
    tot_use1.loc[tot_use1.use_type == 'stockwater', 'use_type_reg'] = 'irrigation'
    tot_use1.loc[tot_use1.use_type == 'water_supply', 'use_type_reg'] = 'other'
    tot_use1.loc[tot_use1.use_type == 'industry', 'use_type_reg'] = 'other'
    tot_use1.loc[tot_use1.use_type == 'hydroelectric', 'use_type_reg'] = 'other'

    ### Create blank multiindex dataframe for both the monthly data and the regressions
    levels0 = et_mon_sel.catch_grp.unique()
    levels1 = ['allo_vol.day', 'usage_vol.day', 'usage_allo_ratio', 'tot_usage_vol.day', 'usage_rate_m3.s']
    levels2 = ['irrigation', 'other', 'tot']
    labels0 = repeat(range(len(levels0)), 8).tolist()
    labels1 = [0,0,1,1,2,2,3,4] * len(levels0)
    labels2 = [0,1,0,1,0,1,2,2] * len(levels0)
    col_index = MultiIndex(levels=[levels0, levels1, levels2], labels=[labels0, labels1, labels2])

    nat_mon = DataFrame(nan, columns=col_index, index=days)

    reg_levels1 = ['irrigation', 'other']
    reg_cols = ["Slope", "Intercept", "R2", "NRMSE", "p-value"]
    reg_row_index = MultiIndex.from_product([levels0, reg_levels1])

    reg_df = DataFrame(nan, columns=reg_cols, index=reg_row_index)

    ### Loop through all sites/files
    usage_et_est =[]
    grp1 = et_mon_sel_rec.groupby('catch_grp')
    for fi, catch_grp in grp1:

        ### Select necessary data
#        usage2 = tot_use1[tot_use1.catch_grp == fi]
#        et_mon_sel2 = et_mon_sel[fi]
#        et_mon_sel_rec2 = et_mon_sel_rec[fi]
        et_mon_sel_rec2 = catch_grp.drop('catch_grp', axis=1).set_index('time')
        tot_use2 = tot_use1.loc[tot_use1.catch_grp == fi, ['crc', 'wap', 'time', 'use_type', 'allo_block', 'take_type', 'use_type_reg', 'sd1_150', mon_allo_name, mon_restr_allo_name, ann_restr_allo_name, usage_name]]
#        tot_use2[date_col] = to_datetime(tot_use2[date_col])
#        usage2[date_col] = to_datetime(usage2[date_col])

#        first_date = et_mon_sel2.index[0]
#        last_date = tot_use2.dates.sort_values().unique()[-1]
#        tot_use2 = tot_use2[(tot_use2[date_col] <= last_date) & (tot_use2[date_col] >= first_date)]

        ### Aggregate all WAP usage by month and use type
        tot_use3 = tot_use2[['time', 'use_type_reg', ann_restr_allo_name, usage_name]].groupby(['time', 'use_type_reg']).sum().round().reset_index()
        tot_use3.columns = ['time', 'use_type_reg', 'allo_vol.day', 'usage_vol.day']
        tot_use3 = tot_use3.pivot(index='time', columns='use_type_reg')
        tot_use3 = tot_use3[in1d(tot_use3.index, days)]

        ### Estimate the number of days per each month
        count1 = tot_use3.index.days_in_month

        ### Convert to volumes per day (to remove the differences in the number of days per month)
        tot_use3 = tot_use3.div(count1.values, axis=0).round(2)
        tot_use3.loc[:, tot_use3.sum() == 0] = nan

        ### Calc ratios and rates
        nat_mon[fi].loc[tot_use3.index, tot_use3.columns] = tot_use3.values
        nat_mon[fi]['usage_allo_ratio'] = (nat_mon[fi]['usage_vol.day']/nat_mon[fi]['allo_vol.day']).round(5)
        nat_mon[fi].loc[:, ('tot_usage_vol.day', 'tot')] = nat_mon[fi]['usage_vol.day'].sum(axis=1)
        nat_mon[fi].loc[:, ('usage_rate_m3.s', 'tot')] = (nat_mon[fi]['tot_usage_vol.day', 'tot'] /24/60/60).round(4)

        ### Lin reg to ratios
        sel1 = concat([nat_mon[fi]['usage_allo_ratio'], et_mon_sel_rec2], axis=1)
        sel1.columns = ['irrigation', 'other', et_col]

        reg1 = lin_reg(sel1[et_col], sel1[['irrigation', 'other']])[0][reg_cols]

        ## Put in mean values if missing
        if reg1.loc['other', :].isnull()[0]:
            reg1.loc['other', ['Slope', 'Intercept']] = [0.000241, 0.023962]
        if reg1.loc['irrigation', :].isnull()[0]:
            reg1.loc['irrigation', ['Slope', 'Intercept']] = [0.000452, -0.004869]

        ## Append to main object
        reg_df.loc[fi] = reg1.values

        ### Estimate the historic usage from ET and allocation
        ## Select and prepare data
#        et_mon_site = et_mon[fi]
        et_mon_site = et_mon_sel_rec2.copy()

        ## calcs
        sd_ratio = DataFrame(nan, index=et_mon_site.index, columns=['irrigation', 'other'])
        sd_ratio['irrigation'] = et_mon_site * reg1['Slope'].values[0] + reg1['Intercept'].values[0]
        sd_ratio['other'] = et_mon_site * reg1['Slope'].values[1] + reg1['Intercept'].values[1]
        sd_ratio2 = sd_ratio.reset_index()
        sd_ratio2.columns = ['time', 'irrigation', 'other']

        # Irrigation
        irr_use1 = tot_use2[tot_use2.use_type_reg == 'irrigation']
        irr_use2 = merge(irr_use1, sd_ratio2[['time', 'irrigation']], on='time', how='left')
        irr_use2['usage.mon'] = irr_use2[ann_restr_allo_name] * irr_use2['irrigation']

        # Other
        oth_use1 = tot_use2[tot_use2.use_type_reg == 'other']
        oth_use2 = merge(oth_use1, sd_ratio2[['time', 'other']], on='time', how='left')
        oth_use2['usage.mon'] = oth_use2[ann_restr_allo_name] * oth_use2['other']

        ## Combine
        tot_use4 = concat([irr_use2.drop('irrigation', axis=1), oth_use2.drop('other', axis=1)])

        ## append the recent better sd data
#        tot_use4 = merge(tot_use3, usage2[['crc', 'wap', date_col, 'use_type', 'sd_usage_est']], on=['crc', 'wap', date_col, 'use_type'], how='left')
        tot_use4.loc[tot_use4[usage_name].isnull(), usage_name] = tot_use4.loc[tot_use4[usage_name].isnull(), 'usage.mon'].round(2)
        tot_use5 = tot_use4.drop(['usage.mon', 'use_type_reg'], axis=1)

        ### Combine into one large object
        usage_et_est.append(tot_use5)

    ### Make the final object
    sd_et_est0 = concat(usage_et_est)[['crc', 'take_type', 'allo_block', 'use_type', 'wap', 'time', 'sd1_150', mon_allo_name, mon_restr_allo_name, ann_restr_allo_name, usage_name]]

    ### Put in actual usage
    sd_et_est1 = merge(sd_et_est0, usage_est[['crc', 'take_type', 'allo_block', 'wap', 'time', raw_usage_name, ann_allo_name]], on=['crc', 'take_type', 'allo_block', 'wap', 'time'], how='left')
    usage_index = sd_et_est1[raw_usage_name].notnull()
    sd_et_est1.loc[usage_index, usage_name] = sd_et_est1.loc[usage_index, raw_usage_name]

    ### Calc SD usage
    sd_et_est1[sd_usage_name] = (sd_et_est1['sd1_150'] * 0.01 * sd_et_est1[usage_name]).round(2)
    sd_et_est1.loc[(sd_et_est1.take_type == 'Take Surface Water'), sd_usage_name] = (sd_et_est1.loc[(sd_et_est1.take_type == 'Take Surface Water'), usage_name]).round(2)

    ## Remove negtive values and others...
    sd_et_est1[sd_et_est1[usage_name] < 0] = 0
    sd_et_est1 = sd_et_est1[sd_et_est1.crc != 0]
    sd_et_est1 = sd_et_est1.drop('sd1_150', axis=1)

    ### Export data and return
    if export:
        save_df(nat_mon, export_mon_path, index=True)
        save_df(reg_df, export_reg_path, index=True)
        save_df(sd_et_est1, export_sd_est_path, index=False)
    return(sd_et_est1, nat_mon, reg_df)


def w_use_proc(ht_use_hdf=r'S:\Surface Water\shared\base_data\usage\ht_usage_daily.h5', export=False, export_path='usage_daily_all.h5'):
    """
    Function to process the water use data.
    """
    from core.ecan_io import rd_sql
    from pandas import to_datetime, read_hdf
    from core.ecan_io.SQL_databases import sql_arg

    #### Import data
    sql1 = sql_arg()
    wus = rd_sql(**sql1.get_dict('wus_day'))
    wus.loc[:, 'usage'] = wus.loc[:, 'usage'].round(2)
    wus.loc[:, 'wap'] = wus.loc[:, 'wap'].str.upper().str.replace(',', '')
    wus.loc[:, 'date'] = to_datetime(wus.loc[:, 'date'])
    wus1 = wus.set_index(['wap', 'date'])

    ht_use = read_hdf(ht_use_hdf).reset_index()
    ht_use.columns = ['wap', 'date', 'usage']
    ht_use.loc[:, 'wap'] = ht_use.loc[:, 'wap'].str.upper().str.replace(',', '')
    ht_use.loc[:, 'date'] = to_datetime(ht_use.loc[:, 'date'])
    ht_use1 = ht_use.set_index(['wap', 'date'])

    #### Merge WUS with new data
    use_daily = wus1.combine_first(ht_use1)

    #### Aggregate waps together
    use_daily2 = use_daily.groupby(level=['wap', 'date']).sum().usage

    #### Export data
    if export:
        use_daily2.to_hdf(export_path, key='daily_usage', mode='w')
    return(use_daily2)

### Errors


def allo_errors(takes, wap, dates, zone, zone_add, takes_names, wap_names, dates_names, zone_names, zone_add_names, irr_names, stock_names, ind_names, pub_names, irr_par=None, irr_par_names=None, export_path=''):
    """
    Function to determine the errors in the database.
    """
    from os import path
    from numpy import in1d
    from pandas import merge, to_datetime, Timestamp
    from core.allo_use import wqn10

    ### Rename columns

    takes.columns = takes_names
    wap.columns = wap_names
    dates.columns = dates_names
    zone.columns = zone_names
    zone_add.columns = zone_add_names
    if irr_par is not None:
        irr_par.columns = irr_par_names

    ### Reorganize and clean the data

    ## Make copies
    takes1a = takes.copy()
    wap1a = wap.copy()
    dates1a = dates.copy()
    zone1a = zone.copy()
    zone_add1a = zone_add.copy()

    ## Select SW and GW takes
    take_types = ['Take Groundwater', 'Take Surface Water']
    takes1b = takes1a[in1d(takes1a.take_type, take_types)]
    wap2 = wap1a[in1d(wap1a.take_type, take_types)]

    ## merge necessary tables

    takes1c = merge(takes1b, dates1a, how='left', on='crc')
    takes2 = merge(takes1c, wap2, how='left', on=['crc', 'take_type'])

    ## select records without max rates and zero values

    err_max_rate_allo = takes1c[(takes1c.max_rate <= 0) | takes1c.max_rate.isnull()]
    err_max_rate_wap = wap2[(wap2.max_rate_wap <= 0) | wap2.max_rate_wap.isnull()]

    err_ret_period_allo = takes1c[(takes1c.max_vol > 0) & ((takes1c.return_period <= 0) | (takes1c.return_period.isnull()))]

    ## Select records that have problems with the 'from' and 'to' dates

    nan_from_to = takes1c[takes1c.from_date.isnull() | takes1c.to_date.isnull()]
    no_wap = takes1c[~in1d(takes1c.crc, wap2.crc)]

    ## Export data
    err_max_rate_allo.to_csv(path.join(export_path, 'err_max_rate_allo.csv'), index=False)
    err_max_rate_wap.to_csv(path.join(export_path, 'err_max_rate_wap.csv'), index=False)
    err_ret_period_allo.to_csv(path.join(export_path, 'err_ret_period_allo.csv'), index=False)
    nan_from_to.to_csv(path.join(export_path, 'err_from_to.csv'), index=False)
    no_wap.to_csv(path.join(export_path, 'err_no_wap.csv'), index=False)

    return()



####################################################
#### Archive


#def allo_use_proc(allo_ts_mon, usage, crc_wap, export=False, export_path='allo_use_ts_mon.csv'):
#    """
#    Function to merge the usage data with the allocation data.
#    """
#    from pandas import merge, to_datetime, DateOffset
#    from core.ecan_io import rd_sql
#    from core.ts import grp_ts_agg
#
#    #### Load data and reformat
##    use1 = rd_sql(code='wus_mon')
##    date_mon = to_datetime(use1.year.astype(str) + '-' + use1.month.astype(str)) + DateOffset(months=1) - DateOffset(days=1)
##    use1['dates'] = date_mon
##    use1 = use1.drop(['year', 'month'], axis=1)
#
#    #### Process the usage data
#    use1 = grp_ts_agg(usage, 'wap', 'date', 'M', 'sum')
##    crc_use1 = use1[use1.wap.str.contains('CRC')]
#    use1.columns = ['wap', 'dates', 'usage']
#    use1.loc[:, 'usage'] = use1.loc[:, 'usage'].round(2)
#    use1.loc[:, 'dates'] = use1.loc[:, 'dates'].astype(str)
#
#    #### Process allo data for the distribution of usage data
#    allo_ts = allo_ts_mon.copy()
#    crc_wap.loc[:, 'wap'] = crc_wap.loc[:, 'wap'].str.upper()
#
#    ### Use type
#    allo_ts['use_type_sum']  = allo_ts.groupby(['crc', 'dates', 'take_type'])['mon_vol'].transform(sum)
#    allo_ts['use_type_ratio'] = allo_ts.mon_vol / allo_ts.use_type_sum
#
#    #### Combine the allo data with the crc_wap
#    crc_wap_dates = merge(allo_ts_mon[['crc', 'dates', 'take_type']], crc_wap, on=['crc', 'take_type']).drop_duplicates()
#    allo_wap1 = merge(crc_wap_dates, allo_ts, on=['crc', 'dates', 'take_type'])
#
#    ### multiple consents on the same month with the same WAP
#    allo_wap1['dup_crc_sum']  = allo_wap1.groupby(['wap', 'dates', 'take_type'])['use_type_sum'].transform(sum)
#    allo_wap1['dup_crc_ratio'] = allo_wap1.use_type_sum / allo_wap1.dup_crc_sum
#
#    ### Combine use with allo/wap
#    allo_use1 = merge(allo_wap1, use1, on=['wap', 'dates'])
#    allo_use1['usage_split'] = allo_use1.dup_crc_ratio * allo_use1.use_type_ratio * allo_use1.usage
#    allo_use2 = allo_use1[['crc', 'dates', 'take_type', 'use_type', 'usage_split']]
#
#    #### Aggregate to crc, month, take_type, and use_type
#    allo_use3 = allo_use2.groupby(['crc', 'dates', 'take_type', 'use_type']).sum().reset_index()
#
#    #### Bring back the allocation
#    allo_use4 = merge(allo_ts_mon, allo_use3, on=['crc', 'dates', 'take_type', 'use_type'], how='right')
#
#    #### Remove oddities
#    allo_use4 = allo_use4[allo_use4.usage_split.notnull()]
#
#    #### Export and return
#    allo_use4.rename(columns={'usage_split': 'usage'}, inplace=True)
#
#    if export:
#        allo_use4.to_csv(export_path, index=False)
#    return(allo_use4)



#def allo_loc_proc(allo, export_shp='crc_loc.shp', export_waps='crc_wap.csv'):
#    """
#    Function to create a shapefile of consent locations. These locations are not the WAPs! They are the locations of the consent owner.
#    """
#    from numpy import in1d, append, nan
#    from pandas import merge, concat
#    from core.ecan_io import rd_sql
#    from core.spatial import xy_to_gpd
#
#    #### Load data
#    crc_loc = rd_sql(code='crc_details_acc_gis')[['crc', 'geometry']]
#
#    waps = rd_sql(code='waps_acc').dropna().drop_duplicates()
#
#    waps_new = rd_sql(code='crc_wap_act_acc')[['crc', 'wap', 'max_rate_wap', 'take_type']].dropna().drop_duplicates()
#
#    crc_rel = rd_sql(code='crc_relation')
#
#    wap_loc = rd_sql(code='waps_details')[['wap', 'NZTMX', 'NZTMY']]
#
#    sd1 = rd_sql(code='sd')[['wap', 'sd1_150']]
#
#    #### Modify the time columns to be strings
#    allo.loc[:, 'from_date'] = allo.loc[:, 'from_date'].astype(str)
#    allo.loc[:, 'to_date'] = allo.loc[:, 'to_date'].astype(str)
#
#    ### Extract all crc
#    crc_set1 = allo.crc.unique()
#
#    #### WAPs
#    ### Combine both waps tables and take the max WAP
#    waps_all1 = merge(waps, waps_new, on=['crc', 'wap', 'take_type'], how='outer')
#    waps_all1.loc[:, 'wap'] = waps_all1.loc[:, 'wap'].str.upper()
#    waps_all1.loc[:, 'wap'] = waps_all1.loc[:, 'wap'].str.replace(' ', '')
#    waps_all1['max_rate'] = waps_all1.max_rate_wap_x
#    null_index = waps_all1['max_rate'].isnull() | (waps_all1['max_rate'] == 0)
#    waps_all1.loc[null_index, 'max_rate'] = waps_all1.loc[null_index, 'max_rate_wap_y']
#    waps_all2 = waps_all1[['crc', 'wap', 'max_rate', 'take_type']].groupby(['crc']).max().reset_index()
#    waps_all2 = waps_all2[in1d(waps_all2.crc, crc_set1)]
#
#    ### Determine missing WAPs
#    crc_mis1 = crc_set1[~in1d(crc_set1, waps_all2.crc.unique())]
#
#    wap_est = []
#    for i in crc_mis1:
#        crc1 = [i]
#        val1 = []
#        while (len(crc1) > 0) & (len(val1) == 0):
#            val1 = waps_all2[in1d(waps_all2.crc, crc1)]
#            crc1 = crc_rel[in1d(crc_rel.crc, crc1)].crc_child.values
#        val1['crc'] = i
#        wap_est.append(val1)
#    wap_est1 = concat(wap_est)
#
#    ## Combine the missing WAPs with the main WAPs object
#    waps_all3 = concat([waps_all2, wap_est1]).drop_duplicates()
#    crc_mis2 = waps_all3[~in1d(waps_all3.wap, wap_loc.wap)].crc.unique()
#
#    waps_all = concat([waps_all1, wap_est1]).drop_duplicates()[['crc', 'wap', 'take_type']]
#    waps_all.loc[:, 'wap'] = waps_all.loc[:, 'wap'].str.upper()
#    waps_all.loc[:, 'wap'] = waps_all.loc[:, 'wap'].str.replace(' ', '')
#    waps_all3 = merge(waps_all3, sd1, on='wap', how='left')
#    waps_all3.loc[waps_all3.sd1_150 == 0, 'sd1_150'] = nan
#
#    ## Second run of missing WAPs
#    crc_mis3 = crc_set1[~in1d(crc_set1, waps_all3.crc.unique())]
#    crc_mis_all = append(crc_mis2, crc_mis3)
#    allo_mis1 = crc_loc[in1d(crc_loc.crc, crc_mis_all)]
#
#    ### Convert x and y to geometry and combine all crc locations
#    waps_all_loc1 = merge(waps_all3, wap_loc, on='wap')
#    allo_loc1 = xy_to_gpd(waps_all_loc1, ['crc', 'sd1_150'], 'NZTMX', 'NZTMY')
#    allo_loc_all = concat([allo_loc1, allo_mis1])
#
#    #### Make final object with stats and locations
#    allo_xy = merge(allo_loc_all, allo, on='crc').drop_duplicates(subset=['crc', 'take_type', 'use_type'])
#
#    #### Save data
#    allo_xy.to_file(export_shp)
#    waps_all.to_csv(export_waps, index=False)
#    return([allo_xy, waps_all])


#
#def allo_ts_proc(allo, end_date='2016-06-30', export=True, export_ann_path='allo_use_yr1.csv', export_mon_path='allo_use_mon1.csv'):
#    from pandas import concat, merge, to_datetime
#    from core.ts import grp_ts_agg
#    """
#    Function set to take the processed allocation data and convert it to a monthly time series. Then the monthly time series is combined with the water usage data and summarized by water year.
#    """
#
#
#    def allo_yr_mon(df, from_col='from_date', to_col='to_date', end_date='2016-06-30'):
#        """
#        Function that converts the allocation data to a monthly time series.
#        """
#        from pandas import DataFrame, DateOffset, date_range, Timestamp
#
#        end_date = Timestamp(end_date)
#        base_ratio = 1.0/12
#        lst1 =[]
#        for i in df.index:
#            date1 = df.loc[i, from_col]
#            date2 = df.loc[i, to_col]
#            day1 = date1.day
#            day2 = date2.day
#            if date1 <= end_date:
#                if date2 > end_date:
#                    date2 = end_date
#                    last1 = 1
#                else:
#                    last1 = float(day2)/date2.days_in_month
#                    date2 = date2 + DateOffset(months=1) - DateOffset(days=1)
#
#                dates1 = date_range(date1, date2, freq='M')
#                first1 = (date1.days_in_month - float(day1))/date1.days_in_month
#
#                if len(dates1) > 0:
#                    df2 = DataFrame({'dates': dates1, 'mon_ratio': [0] * len(dates1)})
#                    df2.loc[:, 'mon_ratio'] = round(base_ratio, 3)
#                    df2.loc[0, 'mon_ratio'] = round(first1 * df2.loc[0, 'mon_ratio'], 3)
#                    df2.loc[len(df2) - 1, 'mon_ratio'] = round(last1 * df2.loc[len(df2) - 1, 'mon_ratio'], 3)
#                    df2.loc[:, 'crc'] = df.loc[i, 'crc']
#                    df2.loc[:, 'take_type'] = df.loc[i, 'take_type']
#                    df2.loc[:, 'use_type'] = df.loc[i, 'use_type']
#
#                    lst1.append(df2)
#
#        df3 = concat(lst1)
#
#        return(df3)
#
#
##    def w_resample(df):
##        df.index = df.dates
##        df_grp = df.resample('A-JUN')
##        df2 = df_grp[['mon_vol', 'mon_vol_wqn']].sum().round(2)
##        return(df2)
#
#    #### Make time series out of consent dates
#    allo.loc[:, 'from_date'] = to_datetime(allo.from_date, errors='coerce')
#    allo.loc[:, 'to_date'] = to_datetime(allo.to_date, errors='coerce')
#    allo_ts1 = allo_yr_mon(allo, end_date=end_date)
#
#    #### Append annual allo volumes and calc monthly volumes
#    allo_ts2 = merge(allo_ts1, allo[['crc', 'take_type', 'use_type', 'ann_vol', 'ann_vol_wqn']], on=['crc', 'take_type', 'use_type'])
#    allo_ts2['mon_vol'] = allo_ts2['mon_ratio'] * allo_ts2['ann_vol']
#    allo_ts2['mon_vol_wqn'] = allo_ts2['mon_ratio'] * allo_ts2['ann_vol_wqn']
#
#    #### drop unneeded columns
#    allo_ts3 = allo_ts2.drop(['mon_ratio', 'ann_vol', 'ann_vol_wqn'], axis=1)
#
#    #### Aggregate to yearly
#    allo_ts4 = grp_ts_agg(allo_ts3, ['crc', 'take_type', 'use_type'], 'dates', 'A-JUN', 'sum').reset_index()
#
#    #### Save data
#    if export:
#        allo_ts3.to_csv(export_mon_path, index=False)
#        allo_ts4.to_csv(export_ann_path, index=False)
#
#    return([allo_ts4, allo_ts3])


def hist_sd_use2(usage_est, allo_gis, vcn_grid2, vcn_data_path, date_col='date', export=True, export_mon_path='sd_mon_vol.csv', export_reg_path='sd_reg.csv', export_sd_est_path='sd_est_all_mon_vol.csv'):
    """
    Function to create stream depletion for all consents from monthly sums of usage and ET.
    """
    from pandas import read_csv, to_datetime, datetime, merge, MultiIndex, DataFrame, date_range, IndexSlice, concat
    from numpy import in1d, array, repeat, nan
    from os import listdir, path
    from fnmatch import filter
    from re import findall
    from core.stats import lin_reg
    from core.misc import rd_dir
    from core.ecan_io import rd_vcn
    from core.ts import w_resample

    ### Name parameters
    usage_name = 'usage_est'
    sd_usage_name = 'sd_usage'

    ### Read in data
    vcn_sites = vcn_grid2[vcn_grid2.ecan_id.notnull()]
    vcn_sites_lst = vcn_sites.ecan_id.astype('int32').values
    vcn_data = rd_vcn(data_dir=vcn_data_path, select=vcn_sites_lst, data_type='ET')
    vcn_sites2 = vcn_sites[in1d(vcn_sites.ecan_id, vcn_data.columns)].sort_values('ecan_id')

    ### Resample the ET series to monthly sums
    vcn_data_agg = vcn_data.groupby(vcn_sites2.catch_grp.values, axis=1).mean()
    et_mon = w_resample(vcn_data_agg, period='month', export=False)
    et_mon = et_mon.dropna(how='all')

    ### Filter only the SW sites with SD and surface water sites
    usage1 = usage_est[(usage_est.sd1_150.notnull()) | (usage_est.take_type == 'Take Surface Water')]
    allo_gis2 = allo_gis[allo_gis.sd1_150.notnull() | (allo_gis.take_type == 'Take Surface Water')].drop('geometry', axis=1).drop_duplicates(subset=['crc', 'take_type', 'use_type'])

    first_date = et_mon.index[0]
    last_date = usage1[usage1.usage_est.notnull()].dates.sort_values().values[-1]
    usage1.loc[:, date_col] = to_datetime(usage1[date_col])
    usage1 = usage1[(usage1[date_col] <= last_date) & (usage1[date_col] >= first_date)]

    ### Estimate the number of days per each month
    days1 = to_datetime(usage1.dates.unique()).sort_values()
    days = days1[in1d(days1, et_mon.index)]
#    count1 = days.days_in_month

    ### Merge allo/use with the spatial info
    tot_use1 = merge(usage1, allo_gis2, how='left')
    tot_use1 = tot_use1[tot_use1.catch_grp.notnull()]

    ### Calc SD usage
#    tot_use1['sd_usage'] = tot_use1['sd1_150'] * 0.01 * tot_use1['usage_est']
#    tot_use1.loc[(tot_use1.take_type == 'Take Surface Water'), 'sd_usage'] = tot_use1.loc[(tot_use1.take_type == 'Take Surface Water'), 'usage_est']

    ### Select ET data that coincide with the usage data
    et_mon_sel = et_mon.loc[:, in1d(et_mon.columns, tot_use1.catch_grp.unique()) & (et_mon.columns != '')]
    et_mon_sel_rec = et_mon_sel.loc[days, :]

    ### Re-group the use types
    tot_use1['use_type_reg'] = tot_use1['use_type']
    tot_use1.loc[tot_use1.use_type == 'stockwater', 'use_type_reg'] = 'irrigation'
    tot_use1.loc[tot_use1.use_type == 'water_supply', 'use_type_reg'] = 'other'
    tot_use1.loc[tot_use1.use_type == 'industry', 'use_type_reg'] = 'other'
    tot_use1.loc[tot_use1.use_type == 'hydroelectric', 'use_type_reg'] = 'other'

    ### Create blank multiindex dataframe for both the monthly data and the regressions
    levels0 = et_mon_sel.columns.tolist()
    levels1 = ['allo_vol.day', 'usage_vol.day', 'usage_allo_ratio', 'tot_usage_vol.day', 'usage_rate_m3.s']
    levels2 = ['irrigation', 'other', 'tot']
    labels0 = repeat(range(len(levels0)), 8).tolist()
    labels1 = [0,0,1,1,2,2,3,4] * len(levels0)
    labels2 = [0,1,0,1,0,1,2,2] * len(levels0)
    col_index = MultiIndex(levels=[levels0, levels1, levels2], labels=[labels0, labels1, labels2])

    nat_mon = DataFrame(nan, columns=col_index, index=days)

    reg_levels1 = ['irrigation', 'other']
    reg_cols = ["Slope", "Intercept", "R2", "NRMSE", "p-value"]
    reg_row_index = MultiIndex.from_product([levels0, reg_levels1])

    reg_df = DataFrame(nan, columns=reg_cols, index=reg_row_index)

    ### Loop through all sites/files
    usage_et_est =[]
    for fi in et_mon_sel.columns.get_level_values(0):

        ### Select necessary data
#        usage2 = tot_use1[tot_use1.catch_grp == fi]
#        et_mon_sel2 = et_mon_sel[fi]
        et_mon_sel_rec2 = et_mon_sel_rec[fi]
        tot_use2 = tot_use1.loc[tot_use1.catch_grp == fi, ['crc', date_col, 'use_type', 'take_type', 'use_type_reg', 'sd1_150', 'mon_vol', 'up_allo_m3', 'ann_up_allo', usage_name]]
        tot_use2[date_col] = to_datetime(tot_use2[date_col])
#        usage2[date_col] = to_datetime(usage2[date_col])

#        first_date = et_mon_sel2.index[0]
#        last_date = tot_use2.dates.sort_values().unique()[-1]
#        tot_use2 = tot_use2[(tot_use2[date_col] <= last_date) & (tot_use2[date_col] >= first_date)]

        ### Aggregate all WAP usage by month and use type
        tot_use3 = tot_use2[[date_col, 'use_type_reg', 'ann_up_allo', usage_name]].groupby([date_col, 'use_type_reg']).sum().round(2).reset_index()
        tot_use3.columns = [date_col, 'use_type_reg', 'allo_vol.day', 'usage_vol.day']
        tot_use3 = tot_use3.pivot(index=date_col, columns='use_type_reg')
        tot_use3 = tot_use3[in1d(tot_use3.index, days)]

        ### Estimate the number of days per each month
        count1 = tot_use3.index.days_in_month

        ### Convert to volumes per day (to remove the differences in the number of days per month)
        tot_use3 = tot_use3.div(count1, axis=0).round(2)
        tot_use3.loc[:, tot_use3.sum() == 0] = nan

        ### Calc ratios and rates
        nat_mon[fi].loc[tot_use3.index, tot_use3.columns] = tot_use3.values
        nat_mon[fi]['usage_allo_ratio'] = (nat_mon[fi]['usage_vol.day']/nat_mon[fi]['allo_vol.day']).round(5)
        nat_mon[fi].loc[:, ('tot_usage_vol.day', 'tot')] = nat_mon[fi]['usage_vol.day'].sum(axis=1)
        nat_mon[fi].loc[:, ('usage_rate_m3.s', 'tot')] = (nat_mon[fi]['tot_usage_vol.day', 'tot'] /24/60/60).round(4)

        ### Lin reg to ratios
        sel1 = concat([nat_mon[fi]['usage_allo_ratio'], et_mon_sel_rec2], axis=1)
        sel1.columns = ['irrigation', 'other', 'et']

        reg1 = lin_reg(sel1['et'], sel1[['irrigation', 'other']])[0][reg_cols]

        ## Put in mean values if missing
        if reg1.loc['other', :].isnull()[0]:
            reg1.loc['other', ['Slope', 'Intercept']] = [0.000241, 0.023962]
        if reg1.loc['irrigation', :].isnull()[0]:
            reg1.loc['irrigation', ['Slope', 'Intercept']] = [0.000452, -0.004869]

        ## Append to main object
        reg_df.loc[fi] = reg1.values

        ### Estimate the historic usage from ET and allocation
        ## Select and prepare data
        et_mon_site = et_mon[fi]

        ## calcs
        sd_ratio = DataFrame(nan, index=et_mon_site.index, columns=['irrigation', 'other'])
        sd_ratio['irrigation'] = et_mon_site * reg1['Slope'].values[0] + reg1['Intercept'].values[0]
        sd_ratio['other'] = et_mon_site * reg1['Slope'].values[1] + reg1['Intercept'].values[1]
        sd_ratio2 = sd_ratio.reset_index()
        sd_ratio2.columns = [date_col, 'irrigation', 'other']

        # Irrigation
        irr_use1 = tot_use2[tot_use2.use_type_reg == 'irrigation']
        irr_use2 = merge(irr_use1, sd_ratio2[[date_col, 'irrigation']], on=date_col, how='left')
        irr_use2['usage.mon'] = irr_use2['ann_up_allo'] * irr_use2['irrigation']

        # Other
        oth_use1 = tot_use2[tot_use2.use_type_reg == 'other']
        oth_use2 = merge(oth_use1, sd_ratio2[[date_col, 'other']], on=date_col, how='left')
        oth_use2['usage.mon'] = oth_use2['ann_up_allo'] * oth_use2['other']

        ## Combine
        tot_use4 = concat([irr_use2.drop('irrigation', axis=1), oth_use2.drop('other', axis=1)])

        ## append the recent better sd data
#        tot_use4 = merge(tot_use3, usage2[['crc', 'wap', date_col, 'use_type', 'sd_usage_est']], on=['crc', 'wap', date_col, 'use_type'], how='left')
        tot_use4.loc[tot_use4[usage_name].isnull(), usage_name] = tot_use4.loc[tot_use4[usage_name].isnull(), 'usage.mon'].round(2)
        tot_use5 = tot_use4.drop(['usage.mon', 'use_type_reg'], axis=1)

        ### Combine into one large object
        usage_et_est.append(tot_use5)

    ### Make the final object
    sd_et_est1 = concat(usage_et_est)[['crc', date_col, 'take_type', 'use_type', 'sd1_150', 'mon_vol', 'up_allo_m3', 'ann_up_allo', usage_name]]

    ### Calc SD usage
    sd_et_est1[sd_usage_name] = (sd_et_est1['sd1_150'] * 0.01 * sd_et_est1[usage_name]).round(2)
    sd_et_est1.loc[(sd_et_est1.take_type == 'Take Surface Water'), sd_usage_name] = (sd_et_est1.loc[(sd_et_est1.take_type == 'Take Surface Water'), usage_name]).round(2)


    ## Remove negtive values and others...
    sd_et_est1[sd_et_est1[sd_usage_name] < 0] = 0
    sd_et_est1 = sd_et_est1[sd_et_est1.crc != 0]
    sd_et_est1 = sd_et_est1.drop('sd1_150', axis=1)

    ### Export data and return
    if export:
        nat_mon.to_csv(export_mon_path)
        reg_df.to_csv(export_reg_path)
        sd_et_est1.to_csv(export_sd_est_path, index=False)
    return([sd_et_est1, nat_mon, reg_df])


def w_use_proc_old(export=False, export_path='usage_daily.h5'):
    """
    Function to process the water use data.
    """
    from core.ecan_io import rd_sql
    from pandas import merge, concat, to_datetime

    #### Import data
    wus = rd_sql(code='wus_day')
    wus.loc[:, 'usage'] = wus.loc[:, 'usage'].round(2)
    wus.loc[:, 'date'] = to_datetime(wus.loc[:, 'date'])

    ht_use = rd_sql('SQL2012TEST01', 'WaterTake', 'UsageReading', ['UsageWap', 'Date', 'Value'])
    ht_use.columns = ['wap_id', 'date', 'usage']
    ht_use_id = rd_sql('SQL2012TEST01', 'WaterTake', 'UsageWap', ['Id', 'Name'])
    ht_use_id.columns = ['wap_id', 'wap']

    #### Process WAP/CRC IDs
    ht_use_id.loc[:, 'wap'] = ht_use_id.wap.str.replace('[:\.]', '/')
#    ht_use_id.loc[ht_use_id.Name == 'L35183/580-M1', 'Name'] = 'L35/183/580-M1' What to do with this one?
    ht_use_id.loc[ht_use_id.wap == 'L370557-M1', 'wap'] = 'L37/0557-M1'
    ht_use_id.loc[ht_use_id.wap == 'L370557-M72', 'wap'] = 'L37/0557-M72'
    ht_use_id = ht_use_id[~ht_use_id.wap.str.contains(' ')]
    ht_use_id.loc[:, 'wap'] = ht_use_id.wap.str.split('-', expand=True)[0]
    ht_use_id = ht_use_id[ht_use_id.wap.str.contains('\d\d\d')]
    ht_use_id.loc[:, 'wap'] = ht_use_id.loc[:, 'wap'].str.upper()
    wus.loc[:, 'wap'] = wus.loc[:, 'wap'].str.upper()

    #### Merge ht use with IDs
    ht_use2 = merge(ht_use, ht_use_id, on='wap_id').drop('wap_id', axis=1)

    #### Merge WUS with new data
    use_daily = concat([wus, ht_use2])

    #### Aggregate waps together
    use_daily2 = use_daily.groupby(['wap', 'date']).sum().reset_index()

    #### Export data
    if export:
        use_daily2.to_hdf(export_path, key='daily_usage', mode='w')
    return(use_daily2)








