# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 23/07/2017 1:31 PM
"""

from __future__ import division
from waimak_extended_boundry.m_packages.wel_packages import get_wel_spd
from waimak_extended_boundry.m_packages.rch_packages import _get_rch
from waimak_extended_boundry import smt
import pandas as pd
import numpy as np

if __name__ == '__main__':
    new_no_flow = smt.get_no_flow()
    zones = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp),'ZONE_CODE')
    zones[~new_no_flow[0].astype(bool)] = np.nan
    # waimak = 4, chch_wm = 7, selwyn=8 , chch_wm chch_formation = 9
    w_idx = np.isclose(zones,4)
    c_idx = np.isclose(zones,7)
    s_idx = np.isclose(zones,8)
    all_idx = np.isfinite(zones)

    outdata = pd.DataFrame(columns=['waimak','selwyn','chch_wm', 'total'])
    well_data = get_wel_spd(smt.wel_version).rename(columns={'row':'i','col':'j', 'layer':'k'})
    well_data['bc_type'] = 'well'

    races = smt.df_to_array(well_data.loc[well_data.type=='race'],'flux')
    wells = smt.df_to_array(well_data.loc[well_data.type=='well'],'flux')
    wells_irr = smt.df_to_array(well_data.loc[(well_data.type=='well') & (well_data.use_type=='irrigation-sw')],'flux')
    wells_oth = smt.df_to_array(well_data.loc[(well_data.type=='well') & (well_data.use_type=='other')],'flux')
    rivers = smt.df_to_array(well_data.loc[well_data.type=='river'],'flux')
    bflux = smt.df_to_array(well_data.loc[well_data.type.str.contains('boundry_flux')],'flux')
    rch = _get_rch()*200*200

    for dat, name in zip([races,wells,wells_irr,wells_oth,rivers, bflux,rch],['races','wells','wells_irr','wells_other','wrivers', 'bflux','rch']):
        for idx, zone in zip([w_idx,c_idx,s_idx,all_idx],['waimak','chch_wm', 'selwyn', 'total']):
            outdata.loc[name,zone] = np.nansum(dat[idx])


    print(outdata/86400)



