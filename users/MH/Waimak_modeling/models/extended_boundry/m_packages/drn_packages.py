"""
Author: matth
Date Created: 20/06/2017 11:58 AM
"""

from __future__ import division
from core import env
import flopy
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.model_tools import get_drn_samp_pts_dict, get_base_drn_cells
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from sfr2_packages import _get_reach_data
from wel_packages import get_wel_spd


def create_drn_package(m, wel_version, reach_version):
    drn_data = _get_drn_spd(wel_version=wel_version, reach_v=reach_version).loc[:,
               ['k', 'i', 'j', 'elev', 'cond']].to_records(False)
    flopy.modflow.mfdrn.ModflowDrn(m,
                                   ipakcb=740,
                                   stress_period_data={0: drn_data},
                                   unitnumber=710)


def _get_drn_spd(reach_v, wel_version):  # todo add pickle at some point


    # load original drains
    drn_data = pd.DataFrame(get_base_drn_cells()) #todo the orginal drain elevations are based off of the reseampeld topo DEM, whihc means theya re as much as +- 8m off

    # take away the cust ones (duplication to id from model where)
    drn_dict = get_drn_samp_pts_dict()
    drn_to_remove = pd.DataFrame(smt.model_where(drn_dict['num7drain_swaz']), columns=['i', 'j'])

    drn_data = pd.concat((drn_data, drn_to_remove))

    idx = ~drn_data.duplicated(['i', 'j'], keep=False)

    drn_data = drn_data.loc[idx]
    drn_data['group'] = None

    # define grouping #todo chat with brioch about this/targets... might want to change grouping for current drains
    index_to_pass = drn_data.index
    for key in drn_dict:
        if ('str' in key or 'swaz' not in key or key == 'num7drain_swaz') and key not in ['waimak_drn', '']:
            continue

        cells = smt.model_where(drn_dict[key])
        drn_to_mark = pd.DataFrame(cells, columns=['i', 'j'], index=range(-1, -len(cells) - 1, -1))
        temp = pd.concat((drn_data, drn_to_mark))
        dup = temp.duplicated(['i', 'j'], False)
        drn_data.loc[dup.loc[index_to_pass], 'group'] = key

    drn_data.loc[np.isclose(drn_data.cond, 585.138611), 'group'] = 'cust_carpet'
    drn_data.loc[np.isclose(drn_data.cond, 5241.530762), 'group'] = 'ash_carpet'
    drn_data.loc[np.isclose(drn_data.cond, 20000.000000), 'group'] = 'chch_carpet'

    # todo the below is a shitty catch all... fix at somepoint
    drn_data.loc[drn_data[
                     'group'].isnull(), 'group'] = 'other'  # this catches a few drain cells that we're not super interested in

    # define zone
    drn_data['zone'] = 'n_wai'  # here this represents the old model
    drn_data.loc[np.in1d(drn_data['group'], ['chch_swaz', 'chch_carpet']), 'zone'] = 's_wai'

    # add a carpet drain south of the waimakariri to loosely represent the low land streams
    # only add drains where there are not other model conditions
    drain_to_add = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/s_carpet_drns.shp".format(smt.sdp),'group',alltouched=True)
    index = np.zeros((smt.rows, smt.cols))

    # wel
    temp = smt.df_to_array(get_wel_spd(wel_version), 'layer', True)[0]
    temp += 1
    index += temp

    # current drain
    temp = smt.df_to_array(drn_data, 'k')
    temp += 1
    index += temp

    # sfr2
    temp = smt.df_to_array((pd.DataFrame(_get_reach_data(reach_v), 'k')))
    temp += 1
    index += temp

    # constant_head/inactive
    temp = smt.get_no_flow(0)
    temp[temp < 0] = 0
    index += ~temp.astype(bool)

    index = index.astype(bool)
    drain_to_add[index] = np.nan

    top = smt.calc_elv_db()[0]
    for val, group in zip([1,3,2,5,4], ['chch_carpet', 'up_lincoln','down_lincoln', 'up_selwyn', 'down_selwyn']):
        temp = pd.DataFrame(smt.model_where(np.isclose(drain_to_add, val)), columns=['i', 'j'])
        temp['k'] = 0
        temp['zone'] = 's_wai'
        temp['group'] = group
        temp['cond'] = 20000
        for i in temp.index:
            row, col = temp.loc['i', 'j']
            temp.loc[i, 'elv'] = top[row, col]

        drn_data = pd.concat(drn_data, temp)

    # add the waimakariri drain up above the bridge
    # set a drain at the top of the ground level
    drain_to_add = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/upper_waimak_drn.shp".format(smt.sdp), 'ID',
                                                 alltouched=True)
    temp = pd.DataFrame(smt.model_where(np.isfinite(drain_to_add)), columns=['i', 'j'])
    temp['k'] = 0
    temp['zone'] = 'n_wai'
    temp['group'] = 'up_waimak'
    temp['cond'] = 20000
    for i in temp.index:
        row, col = temp.loc['i', 'j']
        temp.loc[i, 'elv'] = top[row, col]

    drn_data = pd.concat(drn_data, temp)

    # check for null grouping
    # todo remove all drains in no-flow boundries or constant head
    no_flow = smt.get_no_flow()[0]
    no_flow[no_flow<0]=0
    idxs = pd.DataFrame(smt.model_where(~no_flow.astype(bool)), columns=['i','j'])
    temp_df = pd.concat((drn_data, idxs))
    drn_data = temp_df.loc[~drn_data.duplicated(['i','j'],False)]


    if any(pd.isnull(drn_data['group'])):
        raise ValueError('some groups still null')


if __name__ == '__main__':
    test = _get_drn_spd(1,1)
    print('done')
