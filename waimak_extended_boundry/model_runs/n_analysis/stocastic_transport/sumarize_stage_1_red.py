# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/06/2018 9:03 AM
"""

from __future__ import division
import pandas as pd
from gmp_plus_reductions import get_well_nums_for_group

def get_counts(limit, well_ids, well_counts, datakey, all_data):
    out_num = 0
    temp_data = all_data.loc[well_ids,datakey]
    wells_above = list(temp_data[temp_data>limit].index)

    for wid in wells_above:
        out_num += well_counts[wid]

    return out_num


if __name__ == '__main__':
    wdc_wells = [ # wdc_wells
    'wdc_Kairaki',
    'wdc_Oxford Urban',
    'wdc_Ohoka',
    'wdc_Fernside',
    'wdc_Rangiora',
    'wdc_Mandeville',
    'wdc_West Eyreton',
    'wdc_Kaiapoi',
    'wdc_Waikuku',
    'wdc_Poyntzs Road',
    'wdc_Pegasus',
    'wdc_Cust']

    # private wells
    private = ['Fernside',
    'Flaxton',
    'Horellville',
    'Mandeville',
    'North East Eyrewell_shallow',
    'North West Eyrewell_shallow',
    'Rangiora',
    'Swannanoa_shallow',
    'Waikuku',
    'Woodend - Tuahiwi',
    'West Eyreton_shallow',
    'Clarkville',
    'Cust',
    'Eyreton_deep',
    'Eyreton_shallow',
    'North East Eyrewell_deep',
    'North West Eyrewell_deep',
    'Ohoka_deep',
    'Ohoka_shallow',
    'Springbank',
    'Summerhill',
    'Swannanoa_deep',
    'West Eyreton_deep']


    wells = get_well_nums_for_group()

    all_well_counts = {key: len(val) for key, val in wells.items()}

    data = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\waimakariri_zone\corrected_model_data\n_data_waimak_zone_50ths.csv",
                       index_col=0)

    outdata = pd.DataFrame(index=pd.MultiIndex.from_product((['10% reduction', '30% reduction'],['all land', 'dariy only'])),
                           columns=pd.MultiIndex.from_product((['No of WDC wells', 'No of Private wells'],['> MAV','>7.1 mg/l','> 1/2 MAV'])))
    gmp_data = pd.read_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\all_scens\waimakariri_zone\corrected_model_data\n_data_waimak_zone_50ths.csv", index_col=0)
    for lim, limid in zip([11.3, 7.1, 5.65], ['> MAV', '>7.1 mg/l', '> 1/2 MAV']):
        temp = get_counts(lim, private, well_counts=all_well_counts, datakey='gmp', all_data=gmp_data)
        outdata.loc['gmp',('No of Private wells',limid)] = temp
        temp = get_counts(lim, wdc_wells, well_counts=all_well_counts, datakey='gmp', all_data=gmp_data)
        outdata.loc['gmp', ('No of WDC wells', limid)] = temp

    for red, redid in zip(['red10', 'red20', 'red30'],['10% reduction', '20% reduction', '30% reduction']):
        for scen, scen_id in zip(['kgha_5_','dairy_'],['all land', 'dariy only']):
            for lim, limid in zip([11.3, 7.1, 5.65], ['> MAV', '>7.1 mg/l', '> 1/2 MAV']):
                temp = get_counts(lim,private,well_counts=all_well_counts,datakey=scen + red,all_data=data)
                outdata.loc[(redid,scen_id),('No of Private wells',limid)] = temp
                temp = get_counts(lim,wdc_wells,well_counts=all_well_counts,datakey=scen + red,all_data=data)
                outdata.loc[(redid,scen_id),('No of WDC wells',limid)] = temp

    outdata.to_csv(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\waimakariri_zone_well_summary.csv")
    with open(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\zc_n_sols\stage_1_red\waimakariri_zone_well_summary.csv", 'a') as f:
        f.write('total wdc wells, {}\n'.format(sum([all_well_counts[e]for e in wdc_wells])))
        f.write('total private wells, {}\n'.format(sum([all_well_counts[e]for e in private])))