# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 28/08/2017 4:17 PM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.wel_packages import get_wel_spd
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.well_budget import get_well_budget

def get_well_NSMC_base(version=1):
    all_wells =get_wel_spd(version=version,recalc=True)
    print(get_well_budget(all_wells)/86400)
    all_wells.loc[:,'nsmc_type'] = ''
    #pumping wells
    all_wells.loc[(all_wells.type=='well') & (all_wells.cwms=='chch'),'nsmc_type'] = 'pump_c'
    all_wells.loc[(all_wells.type=='well') & (all_wells.cwms=='selwyn'),'nsmc_type'] = 'pump_s'
    all_wells.loc[(all_wells.type=='well') & (all_wells.cwms=='waimak'),'nsmc_type'] = 'pump_w'

    #selwyn_hillfeds
    all_wells.loc[all_wells.type=='river','nsmc_type'] = 'sriv'

    #races
    all_wells.loc[(all_wells.type == 'race')&(all_wells.zone=='n_wai'),'nsmc_type'] = 'n_race'
    all_wells.loc[(all_wells.type == 'race')&(all_wells.zone=='s_wai'),'nsmc_type'] = 's_race'

    #boundary fluxes
    all_wells.loc[all_wells.type == 'llr_boundry_flux','nsmc_type'] = 'llrzf'
    all_wells.loc[all_wells.type == 'ulr_boundry_flux','nsmc_type'] = 'ulrzf'
    all_wells.loc[all_wells.type == 'boundry_flux','nsmc_type'] = 'nbndf'

    all_wells = all_wells.loc[:,['layer','row','col', 'flux','nsmc_type']]
    return all_wells

def get_template_data():
    groups = ['pump_c', 'pump_s', 'pump_w', 'sriv', 'n_race', 's_race', 'nbndf','llrzf', 'ulrzf']
    pram_id = ['${}$'.format(e) for e in groups]
    out_data = pd.DataFrame(index=groups,data=pram_id)
    return out_data

if __name__ == '__main__':
    well_data = get_well_NSMC_base(3)
    well_data.to_csv(r"C:\Users\MattH\Desktop\to_brioch_2017_10_4/well_data.csv")