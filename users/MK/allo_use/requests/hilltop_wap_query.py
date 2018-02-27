# -*- coding: utf-8 -*-
"""
Script to query the usage data from hilltop xml exports.
"""

from core.ecan_io.hilltop import iter_xml_dir, data_check_fun, proc_use_data, parse_ht_xml, all_data_fun, convert_site_names, rd_hilltop_sites
from pandas import concat, read_csv, to_datetime, DataFrame, Series
from os.path import join, basename
from numpy import nan, all, array_equal
from core.misc import rd_dir

#########################################
#### Parameters

#fpath = {'tel': r'C:\ecan\hilltop\xml_test\tel', 'annual': r'C:\ecan\hilltop\xml_test\annual', 'archive': r'C:\ecan\hilltop\xml_test\archive'}
fpath = r'H:\Data\Telemetry'

export_waps = r'E:\ecan\local\Projects\requests\Ilja\2017-07-25\hilltop_telemetry_waps.csv'
export_waps1 = r'E:\ecan\local\Projects\requests\Ilja\2017-07-25\hilltop_telemetry_waps_with_m.csv'

#########################################
#### Read WAPs in directory

files1 = rd_dir(fpath, 'hts')

data_lst = []
for i in files1:
    f1 = rd_hilltop_sites(join(fpath, i))
    f1['hts_file'] = i
    data_lst.append(f1)

results1 = concat(data_lst)

#########################################
#### Aggregate to WAP without M's

results2 = results1.copy()

results2.loc[:, 'site'] = results2.loc[:, 'site'].apply(lambda x: x.split('-')[0])
results3 = results2.drop_duplicates('site')

#########################################
#### Output

results3.to_csv(export_waps, index=False)
results1.to_csv(export_waps1, index=False)





















########################################
#### Run stats iteration

res1 = []
for i in fpath:
    res_file = iter_xml_dir(fpath[i], stats_fun=data_check_fun, with_xml=True, export=True, export_name=i + '.csv')
    res1.append(res_file)

df_res = concat(res1)

### If runs were already performed

res1 = []
for i in fpath:
    res_file = read_csv(join(fpath[i], i + '.csv'))
    res1.append(res_file)
df_res = concat(res1)

df_res.loc[:, 'end_date'] = to_datetime(df_res.loc[:, 'end_date'])
df_res.loc[:, 'start_date'] = to_datetime(df_res.loc[:, 'start_date'])

#####################################
#### Overall checks

len(df_res.site.unique())

### Negatives

neg_thres = 0.01
neg_waps = df_res[df_res['n_neg/n_data'] >= neg_thres]
len(neg_waps)

# Using a threshold of 0.01 as a cut off for WAPs seems appropriate

### Zeros

zero_thres = 0.95
zero_waps = df_res[df_res['n_zero/n_data'] >= zero_thres]
len(zero_waps)

### Duplicate WAPs

grp1 = df_res.groupby('site')
count1 = grp1['n_data'].count()
(count1 > 1).sum()
count1[count1 > 2]


grp2 = df_res.groupby(['site', 'mtype'])
count2 = grp2['n_data'].count()
(count2 > 1).sum()


t2 = df_res[df_res.site == 'BZ19/0098-M1']
t1 = df_res[df_res.site == 'BX22/0006-M1']

def overlap(df):
    start = df['start_date'].values
    end = df['end_date'].values
    count0 = sum([sum(i < end) for i in start])
    olap = (count0 == len(start)*2) | (count0 > 6)
    return(olap)

ol1 = grp1[['start_date', 'end_date']].apply(overlap)

ol2 = ol1[ol1].index

df = df_res[df_res.site == 'N32/0147-M1']

t1 = df_res[df_res.site.isin(ol2)]
ol_file = t1.groupby('file_name')['n_data'].count()
count_file = df_res.groupby('file_name')['n_data'].count()

ol_tot = concat([ol_file, count_file], axis=1)
ol_tot.columns = ['Overlapping_sites', 'Total_sites']

ol_tot.to_csv(r'C:\ecan\hilltop\xml_test\overlapping_sites.csv')
df_res.to_csv(r'C:\ecan\hilltop\xml_test\all_sites.csv', index=False)

#######################################################
#### Filter through and select the useful sites by file


def overlap(df):
    start = df['start_date'].values
    end = df['end_date'].values
    count0 = sum([sum(i < end) for i in start])
    olap = (count0 == len(start)*2) | (count0 > 6)
    return(olap)


def sel_dup_mtype_sites(df):
    if len(df) == 1:
        return(df)
    else:
        sboo = df.start_date.isin([df.start_date.min()])
        eboo = df.end_date.isin([df.end_date.max()])
        one_rules = sboo & eboo
        if all([sboo, eboo]):
            index1 = [False] * len(df)
            index1[0] = True
            return(df.loc[index1])
        elif any(one_rules):
            return(df.loc[one_rules])
        else:
            idmin = df.start_date.idxmin()
            idmax = df.end_date.idxmax()
            return(df.loc[[idmin, idmax]])


### Must have normal WAP names
df = df_res.copy()
df.rename(columns={'site': 'wap'}, inplace=True)
df.loc[:, 'wap_name'] = df.wap.str.replace('[:\.]', '/')
#    df.loc[df.Name == 'L35183/580-M1', 'Name'] = 'L35/183/580-M1' What to do with this one?
df.loc[df.wap_name == 'L370557-M1', 'wap_name'] = 'L37/0557-M1'
df.loc[df.wap_name == 'L370557-M72', 'wap_name'] = 'L37/0557-M72'
df.loc[df.wap_name.str.contains(' '), 'wap_name'] = nan
df.loc[:, 'wap_name'] = df.wap_name.str.split('-', expand=True)[0]
df.loc[~df.wap_name.str.contains('\d\d\d', na=True), 'wap_name'] = nan
#df.loc[df.wap_name.str.contains('-M'), 'wap_name'] = nan
df.loc[:, 'wap_name'] = df.loc[:, 'wap_name'].str.upper()

### No negative values over 0.01
neg_thres = 0.01
df.loc[df['n_neg/n_data'] > neg_thres, 'wap_name'] = nan

### Remove sites with less than 300 data points
df.loc[df['n_data'] <= 300, 'wap_name'] = nan

### Remove sites with dates earlier than 1990
df.loc[df['start_date'] < '1990', 'wap_name'] = nan

### Output these problem sites
bad1 = df[df.wap_name.isnull()]
bad1.to_csv(export_bad1, index=False)

### Deal with duplicates
df2 = df[df.wap_name.notnull()]
grp1 = df2.groupby(['wap', 'mtype'])
grp1 = df2.groupby('wap')

out2 = []
for name, group in grp1:
    out1 = sel_dup_mtype_sites(group)
    out2.append(out1)
df3 = concat(out2)

df3.to_csv(export_usable, index=False, encoding='utf-8')

##################################################
#### Process water use data

export_name = '_daily'

### Read in data
usites = read_csv(export_usable)

### Select only the file_name, ht site name, and mtype
select = usites[['xml', 'wap', 'mtype']]
select.columns = ['file_name', 'site', 'mtype']

### Run through all sites!
annual_daily = iter_xml_dir(fpath['annual'], stats_fun=proc_use_data, with_xml=True, select=select, export=True, export_name='annual' + export_name + '.csv')
tel_daily = iter_xml_dir(fpath['tel'], stats_fun=proc_use_data, with_xml=True, select=select, export=True, export_name='tel' + export_name + '.csv')
archive_daily = iter_xml_dir(fpath['archive'], stats_fun=proc_use_data, with_xml=True, select=select, export=True, export_name='archive' + export_name + '.csv')

### Load them in if already done
tel_daily = read_csv(join(fpath['tel'], 'tel' + export_name + '.csv'))
archive_daily = read_csv(join(fpath['archive'], 'archive' + export_name + '.csv'))
annual_daily = read_csv(join(fpath['annual'], 'annual' + export_name + '.csv'))

### Combine and export
use1 = concat([tel_daily, archive_daily, annual_daily]).reset_index(drop=True)
use1.to_csv(export_daily_raw, index=False)

### Aggregate sites (first with the 'M's) - use last data set when overlapping


def agg_sites_xml(grp):
    """
    Groupby function to combine overlapping data sets from multiple providers.
    """

    if len(grp.xml.unique()) == 1:
        grp2 = grp[['date', 'val']].copy()
        grp2.set_index('date', inplace=True)
        return(grp2.val)
    else:
        print(grp.name)
        grp2 = grp.copy()
        grp2.set_index('date', inplace=True)
        xmls = grp.groupby('xml')['date'].max().sort_values(ascending=False)
        set1 = grp2.loc[grp2.xml == xmls.index[0], 'val'].groupby(level=0).mean()
        set2 = grp2.loc[grp2.xml == xmls.index[1], 'val'].groupby(level=0).mean()
        grp3 = set1.combine_first(set2)
        return(grp3)

use2 = use1.groupby('site').apply(agg_sites_xml)
use2.to_csv(export_daily_meters)

### Aggregate WAPs regardless of 'M's

use3 = use2.reset_index()
site_names = convert_site_names(use3.site)
use3.loc[:, 'site'] = site_names

use4 = use3.groupby(['site', 'date']).sum().reset_index()
use4.to_csv(export_daily_waps, index=False)

## Pivot

use5 = use4.pivot('date', 'site', 'val')
use5.to_csv(export_daily_waps_cols)






###############################################
#### Testing

wap = 'BY20/0130-M1'
wap = 'K37/0729'
df = df2[df2.wap == wap]


a1 = df1[df1.mtype == 'Flow']
a1.time_res.unique().tolist()



df.iloc[1,1] = df.iloc[2,1]
df.iloc[1,2] = df.iloc[0,2]

sboo = df.start_date.isin([df.start_date.min()])
eboo = df.end_date.isin([df.end_date.max()])


sboo & eboo



dff = DataFrame({'A': np.arange(8), 'B': list('aabbbbcc')})
dff['C'] = np.arange(8)
dff.groupby('B').filter(lambda x: len(x['C']) > 2)



df2[df2.wap == 'swap-M1169']

count2 = df3.groupby(['wap'])['n_data'].count().reset_index()
dups = count2[count2 > 1]

df3.loc[df3.wap.isin(dups), ['file_name', 'start_date', 'end_date', 'wap', 'mtype']]

xml = r'C:\ecan\hilltop\xml_test\annual\Anonymous_Flow.xml'
xml = r'C:\ecan\hilltop\xml_test\tel\WaterOutlook.xml'
select1 = select.loc[select.file_name == 'WaterOutlook.xml', ['site', 'mtype']]
output = parse_ht_xml(xml, proc_use_data, select1)

usites[usites.file_name == 'Anonymous_Flow.xml']
s1 = DataFrame([['N34/0342-M1', 'Flow']], columns=['site', 'mtype'])

t1 = parse_ht_xml(xml, proc_use_data, s1)
t2 = parse_ht_xml(xml, all_data_fun, s1)

data = t2.val


wap1 = 'BV25/0003-M1'
t1 = use1[use1.site == wap1]
t1.groupby(['site', 'date']).sum()

xml = r'C:\ecan\hilltop\xml_test\tel\Harvest.xml'

grp = use1.loc[use1.site == 'BV25/0003-M1', ['date', 'val', 'xml']]
grp = use1.loc[use1.site == 'BW22/0047-M1', ['date', 'val', 'xml']]
grp = use1.loc[use1.site == 'P30/0019-M1', ['date', 'val', 'xml']]

usites[usites.wap_name == 'P30/0019']

s1 = DataFrame([['P30/0019-M1', 'Water Meter']], columns=['site', 'mtype'])

t2 = parse_ht_xml(xml, all_data_fun, s1)


use1 = read_csv(export_daily_waps)
wap = 'J36/0016'

use1[use1.site == wap]


n1 = df_res.groupby('site')['n_data'].count()
n1.index[n1 > 2]

grp = use1[use1.site == 'K38/2012-M1']

df_res[df_res.site == 'K38/2012-M1']




































