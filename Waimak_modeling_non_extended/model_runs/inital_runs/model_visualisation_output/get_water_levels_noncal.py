"""
Author: matth
Date Created: 17/05/2017 3:50 PM
"""

from __future__ import division
import env
import flopy
import numpy as np
import users.MH.Waimak_modeling.model_tools as mt
from users.MH.Waimak_modeling.supporting_data_path import sdp
from core.classes.hydro import hydro
import pandas as pd
from core.ecan_io import rd_sql, sql_db


def main():
    out_path = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/"
                       "Well_BCsand Targets/Steady state data")
    hds = flopy.utils.HeadFile("{}/base_model_runs/base_ss_mf/base_SS.hds".format(sdp))
    head = hds.get_data(hds.get_kstpkper()[-1])
    head[np.isclose(head, -999)] = np.nan
    wells = mt.get_all_well_data()
    with open('{}/Non-target_WL_wells.txt'.format(out_path)) as f:
        well_list = f.readlines()
    well_list = [e.strip() for e in well_list]
    out_wells = wells[np.in1d(wells.index, well_list)]
    out_wells['head'] = np.nan
    out_wells['obs_min'] = np.nan
    out_wells['obs_q25'] = np.nan
    out_wells['obs_q50'] = np.nan
    out_wells['obs_q75'] = np.nan
    out_wells['obs_max'] = np.nan
    out_wells['depth'] = np.nan
    gwl2 = hydro().get_data(mtypes='gwl_m', sites=out_wells.index)
    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.set_index('WELL_NO')

    for well in out_wells.index:
        try:
            ref_level = well_details.loc[well, 'REFERENCE_RL']
            out_wells.loc[well,'depth'] = well_details.loc[well,'DEPTH']
        except:
            ref_level = None

        idx = (out_wells.loc[well, 'layer'], out_wells.loc[well, 'row'], out_wells.loc[well, 'col'])
        out_wells.loc[well, 'head'] = head[idx]

        if pd.isnull(ref_level):
            continue

        try:
            out_wells.loc[well, 'obs_min'] = gwl2.data.loc[:, well].min() + ref_level
            out_wells.loc[well, 'obs_q25'] = gwl2.data.loc[:, well].quantile(0.25) + ref_level
            out_wells.loc[well, 'obs_q50'] = gwl2.data.loc[:, well].quantile(0.50) + ref_level
            out_wells.loc[well, 'obs_q75'] = gwl2.data.loc[:, well].quantile(.75) + ref_level
            out_wells.loc[well, 'obs_max'] = gwl2.data.loc[:, well].max() + ref_level
        except Exception as val:
            print(val)

    out_wells.to_csv('{}/heads_from_ss_for_non_target_wells.csv'.format(out_path))


if __name__ == '__main__':
    main()
    out_path = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/"
                       "Well_BCsand Targets/Steady state data")
    d = pd.read_csv('{}/heads_from_ss_for_non_target_wells.csv'.format(out_path))
    d = d[d['depth']>=50]
    m = mt.wraps.mflow.get_base_mf_ss()
    fig,ax,mmap = mt.plt_default_map(m,0)
    scat = ax.scatter(d['x'],d['y'],c=d['head']-d['obs_q50'],cmap='RdBu',vmin=-15,vmax=15,s=50)
    fig.colorbar(scat, ax=ax, extend='max')
    fig.suptitle('Wells >= 50m')
    fig.savefig('{}/non_target_wells.png'.format(out_path))


