# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 3/10/2017 2:14 PM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import pandas as pd
import numpy as np

def get_well_budget(well_data):
    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp),'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = np.nan
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    w_idx = np.isclose(zones,4)
    c_idx = np.isclose(zones,7)
    s_idx = np.isclose(zones,8)
    all_idx = np.isfinite(zones)
    outdata = pd.DataFrame(columns=['waimak','selwyn','chch_wm', 'total'])
    #well_data = well_data.rename(columns={'row':'i','col':'j', 'layer':'k'})
    well_data['bc_type'] = 'well'

    races = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[well_data.type=='race'])),'flux')
    wells = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[well_data.type=='well'])),'flux')
    wells_irr = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[(well_data.type=='well') & (well_data.use_type=='irrigation-sw')])),'flux')
    wells_oth = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[(well_data.type=='well') & (well_data.use_type=='other')])),'flux')
    rivers = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[well_data.type=='river'])),'flux')
    bflux = smt.df_to_array(pd.DataFrame(smt.convert_well_data_to_stresspd(well_data.loc[well_data.type.str.contains('boundry_flux')])),'flux')

    for dat, name in zip([races,wells,wells_irr,wells_oth,rivers, bflux],['races','wells','wells_irr','wells_other','wrivers', 'bflux']):
        for idx, zone in zip([w_idx,c_idx,s_idx,all_idx],['waimak','chch_wm', 'selwyn', 'total']):
            outdata.loc[name,zone] = np.nansum(dat[idx])

    return outdata
