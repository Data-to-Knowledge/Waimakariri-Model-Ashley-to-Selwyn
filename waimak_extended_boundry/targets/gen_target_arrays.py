# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/07/2017 12:55 PM
"""

from __future__ import division
import env
import numpy as np
import pandas as pd
from waimak_extended_boundry.extended_boundry_model_tools import smt, _get_constant_heads
from waimak_extended_boundry.m_packages.drn_packages import _get_drn_spd
import geopandas as gpd


def gen_drn_target_array():
    drn_data = _get_drn_spd(smt.reach_v, smt.wel_version)
    out_dict = {}
    out_array = np.zeros((smt.rows, smt.cols))
    for i, group in enumerate(set(drn_data['target_group']), 1):
        temp = drn_data.loc[drn_data.target_group == group]
        temp_array = smt.df_to_array(temp, 'k')
        out_array[np.isfinite(temp_array)] = i
        out_dict[i] = group

    kspit = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/kspit.shp".format(smt.sdp), 'Id', True)
    out_array[np.isfinite(kspit) & (np.isclose(out_array, 22) | np.isclose(out_array, 3))] = 45
    out_dict[45] = 'd_kspit'

    return out_array.astype(int), out_dict


def gen_sfr_flow_target_array():
    shp_path = '{}/m_ex_bd_inputs/shp/org_str_flow_targets.shp'.format(smt.sdp)
    target_array = smt.shape_file_to_model_array(shp_path, 'RASTERVALU', True)
    target_array[np.isnan(target_array)] = 0
    num_to_name = {1: 'sfo_c_benn',
                   2: 'sfo_c_tip',
                   3: 'sfo_c_pat',
                   4: 'sfo_c_oxf',
                   5: 'sfo_c_swan',
                   6: 'sfo_1drn',
                   7: 'sfo_2drn',
                   8: 'sfo_3drn',
                   9: 'sfo_4drn',
                   10: 'sfo_5drn',
                   11: 'sfo_7drn',
                   12: 'sfo_c_tlks',
                   13: 'sfo_c_skew',
                   14: 'sfo_e_wolf',
                   15: 'sfo_e_poyn',
                   16: 'sfo_e_down',
                   17: 'sfo_e_seyr'}
    return target_array.astype(int), num_to_name


def gen_sfr_full_we_flux_target_array():
    shp_path = '{}/m_ex_bd_inputs/shp/full_w_e_flux_targets.shp'.format(smt.sdp)
    target_array = smt.shape_file_to_model_array(shp_path, 'targ_code', True)
    target_array[np.isnan(target_array)] = 0
    num_to_name = {2: 'sfx_w_all',
                   1: 'sfx_e_all',
                   3: 'sfx_cd_all'}

    return target_array.astype(int), num_to_name


def gen_sfr_flux_target_array():
    shp_path = '{}/m_ex_bd_inputs/shp/org_str_flux_targets.shp'.format(smt.sdp)
    target_array = smt.shape_file_to_model_array(shp_path, 'GRID_CODE', True)
    target_array[np.isnan(target_array)] = 0
    num_to_name = {10: 'sfx_a1_con',
                   11: 'sfx_c1_swa',
                   12: 'sfx_c2_mil',
                   13: 'sfx_e1_stf',
                   14: 'sfx_w1_cou',
                   15: 'sfx_w2_tom',
                   16: 'sfx_w3_ros',
                   17: 'sfx_w4_mcl',
                   18: 'sfx_w5_wat',
                   19: 'sfx_w6_sh1',
                   22: 'sfx_a2_gol',
                   23: 'sfx_a4_sh1',
                   24: 'sfx_a3_tul',
                   115: 'sfx_3drn',
                   116: 'sfx_4drn',
                   117: 'sfx_5drn',
                   118: 'sfx_6drn',
                   119: 'sfx_2adrn',
                   120: 'sfx_7drn',
                   888: 'sfx_a0_gry'}

    return target_array.astype(int), num_to_name


def gen_constant_head_targets():  # watch if we have constant heads in the sw boundary also note that this is 3d
    chbs = _get_constant_heads()
    shp_path = "{}/m_ex_bd_inputs/shp/coastal_head_target_zones.shp".format(smt.sdp)
    zones = smt.shape_file_to_model_array(shp_path, 'Id', alltouched=True)
    zones[np.isnan(chbs[0])] = 0
    zones[np.isnan(zones)] = 0
    inland = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/inland_zones.shp".format(smt.sdp), 'Id',
                                           alltouched=True)
    zones[np.isfinite(inland) & (np.isclose(zones, 0))] = 9
    zones = np.repeat(zones[np.newaxis, :, :], 11, axis=0)
    no_flow = smt.get_no_flow()
    zones[np.isclose(no_flow, 0)] = 0
    zone_data = gpd.read_file(shp_path)
    zone_data = zone_data.set_index('Id')
    zone_data = zone_data.loc[:, 'name'].to_dict()
    zone_data[9] = 'onshore'

    return zones.astype(int), zone_data


def get_target_group_values():
    # Note that below are in m3/s and then converted to m3/day
    target_group_val = {'chb_ash': -2.4,
                        'chb_chch': -0.8,
                        'chb_cust': -0.35,
                        'chb_sely': 'sel_off',
                        'onshore': None,

                        # drains
                        'd_ash_c': None,
                        'd_ash_est': -0.15,
                        'd_ash_s': None,
                        'd_bul_avon': 'chch_str',
                        'd_bul_styx': 'chch_str',
                        'd_cam_s': None,
                        'd_chch_c': 'chch_str',
                        'd_court_s': None,
                        'd_cust_c': None,
                        'd_dlin_c': 'sel_str',
                        'd_dsel_c': 'sel_str',
                        'd_dwaimak': None,
                        'd_ect': None,
                        'd_greig_s': None,
                        'd_kaiapo_s': None,
                        'd_salt_s': None,
                        'd_taran_s': None,
                        'd_ulin_c': 'sel_str',
                        'd_usel_c': 'sel_str',
                        'd_uwaimak': None,
                        'd_waihora': 'sel_off',
                        'd_waikuk_s': None,
                        'd_kspit': 'sel_off',

                        # from previous shapefile of targets
                        'd_cam_mrsh': -0.17,
                        'd_cam_revl': 0.0,
                        'd_cam_yng': -0.33,
                        'd_cour_nrd': -0.37,
                        'd_emd_gard': -0.18,
                        'd_kairaki': -0.09,
                        'd_kuku_leg': -0.45,
                        'd_nbk_mrsh': -0.66,
                        'd_oho_btch': -0.29,
                        'd_oho_jefs': -0.03,
                        'd_oho_kpoi': -0.19,
                        'd_oho_misc': 0.0,
                        'd_oho_mlbk': -0.15,
                        'd_oho_whit': -0.01,
                        'd_salt_fct': -0.14,
                        'd_salt_top': -0.27,
                        'd_sbk_mrsh': -0.17,
                        'd_sil_harp': -0.36,
                        'd_sil_heyw': -0.39,
                        'd_sil_ilnd': -0.89,
                        'd_smiths': -0.12,
                        'd_tar_gre': -0.11,
                        'd_tar_stok': -0.10,

                        # surface water flow # from previous shapefiles of targets
                        'sfo_1drn': None,
                        'sfo_2drn': None,
                        'sfo_3drn': None,
                        'sfo_4drn': None,
                        'sfo_5drn': None,
                        'sfo_7drn': None,
                        'sfo_c_benn': 0.13,
                        'sfo_c_oxf': 0.99,
                        'sfo_c_pat': None,  # was previously 0 but this meant no data
                        'sfo_c_skew': 1.75,
                        'sfo_c_swan': 0.85,
                        'sfo_c_tip': None,  # was previously 0 but this meant no data
                        'sfo_c_tlks': 1.75,
                        'sfo_e_wolf': 0.0,
                        'sfo_e_poyn': 0.0,
                        'sfo_e_down': 0.0,
                        'sfo_e_seyr': 0.0,

                        # surface water flux #from pervious shapefiles of targets
                        'sfx_a0_gry': 0.2,  # this was added on 19/12/2017 as we found the optimised model lost heaps of water
                        'sfx_a1_con': 'mid_ash_g',
                        'sfx_a2_gol': 'mid_ash_g',
                        'sfx_a3_tul': 'mid_ash_g',
                        'sfx_a4_sh1': -0.40,
                        'sfx_c1_swa': None,  # previously but doesn't add up -0.43,
                        'sfx_c2_mil': None,  # previously but doesn't add up -0.42,
                        'sfx_e1_stf': 1.6,
                        'sfx_w1_cou': -0.9,  # Waimakariri Gorge- Courtenay Road
                        'sfx_w2_tom': 1.1,  # Courtenay Road - Halkett
                        'sfx_w3_ros': 2.8,  # Halkett- Weedons Ross Road
                        'sfx_w4_mcl': 4.7,  # Weedons Ross Road - Crossbank
                        'sfx_w5_wat': 3.2,  # Crossbank - Wrights Cut
                        'sfx_w6_sh1': 1.1,  # Wrights Cut - Old State Highway Bridge

                        # below are the old waimak targets kept here for posterity
                        # these targets assume that the model can model underflow in the
                        # waimakariri river, which we think is implausable given the 30m+
                        # layer 0 thickness
                        # 'sfx_w1_cou': 2.8, # old kept for posterity
                        # 'sfx_w2_tom': 1.1, # old kept for posterity
                        # 'sfx_w3_ros': 2.2, # old kept for posterity
                        # 'sfx_w4_mcl': 5.7, # old kept for posterity
                        # 'sfx_w5_wat': -0.1, # old kept for posterity
                        # 'sfx_w6_sh1': -0.5, # old kept for posterity

                        'sfx_3drn': None,
                        'sfx_2adrn': None,
                        'sfx_4drn': None,
                        'sfx_5drn': None,
                        'sfx_6drn': None,
                        'sfx_7drn': None,
                        'sfx_e_all': 1.99,
                        'sfx_w_all': 12,
                        'sfx_cd_all': None,

                        # groups
                        'mid_ash_g': 5.20,
                        'sel_off': -5.9,  # 1.2-17 # to much range- just observe in NSMC
                        'chch_str': -10,  # 7.5 to 12.5
                        'sel_str': -9.8  # no range
                        }
    for key in target_group_val.keys():
        if isinstance(target_group_val[key], str) or target_group_val[key] is None:
            continue
        target_group_val[key] *= 86400  # convert numbers to m3/day

    return target_group_val


def get_vertical_gradient_targets():
    # load in vert targets
    vert_targets = pd.read_excel(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/Vertical gradient targets updated_use.xlsx"),
        sheetname='data_for_python', index_col=0)

    # load in the row, col options from the bulk head targets sheet
    all_wells = pd.read_csv('{}/all_wells_row_col_layer.csv'.format(smt.sdp), index_col=0)
    vert_targets = pd.merge(vert_targets, all_wells, how='left', left_index=True, right_index=True)

    vert_targets.loc['M35/11937', 'GWL_RL'] = vert_targets.loc[['M35/11937', 'M35/10909'], 'GWL_RL'].mean()
    vert_targets = vert_targets.drop(['M35/10909'])  # this and above are in the same layer

    all_targets = pd.read_csv(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_inc_error.csv"),
        index_col=0)
    head_targets = get_head_targets(True)
    idx = vert_targets.index
    # reviesd error for the NSMc as of 24/10/2017 take the total error and scale it between 1 and 10
    vert_targets.loc[idx, 'weight_2'] = head_targets.loc[idx,'weight_2']
    vert_targets.loc[idx, 'weight'] = 1 / all_targets.loc[idx, 'total_error_m']
    vert_targets.loc[idx, 'GWL_RL'] = all_targets.loc[idx, 'h2o_elv_mean']

    all_targets = vert_targets.loc[:,
                  ['NZTM_x', 'NZTM_y', 'layer', 'GWL_RL', 'weight', 'weight_2', 'row', 'col', 'group', 'Location', 'code']].rename(
        columns={'NZTM_x': 'x', 'NZTM_y': 'y', 'GWL_RL': 'obs', 'row': 'i', 'col': 'j'})
    # return a dataframe: lat, lon, layer, obs, weight?, i,j
    # pull weight from uncertainty
    all_targets = all_targets.dropna()
    all_targets = all_targets.loc[all_targets.duplicated('group', keep=False)]

    return all_targets


def get_head_targets(return_missing = False):
    """

    :param return_missing: boolean if true return 'BW23/0165' (for incorporation of new target weights)
    :return:
    """
    # lat, lon, layer, obs, weigth? i, j
    all_targets = pd.read_csv(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/targets/head_targets/head_targets_2008_inc_error.csv"),
        index_col=0)
    all_targets = all_targets.loc[(all_targets.h2o_elv_mean.notnull()) & (all_targets.row.notnull()) &
                                  (all_targets.col.notnull()) & (all_targets.layer.notnull())]

    # pull out targets for each layer
    min_readings = {0: 5,
                    1: 2,
                    2: 1,
                    3: 1,
                    4: 1,
                    5: 1,
                    6: 1,
                    7: 1,
                    8: 1,
                    9: 1}
    all_targets.loc[:, 'weight'] = 1 / all_targets.loc[:, 'total_error_m']

    outdata = pd.DataFrame()
    for layer in range(smt.layers - 1):  # pull out targets for layers 0-9 layer 10 has no targets
        idx = (all_targets.layer == layer) & (all_targets.readings_nondry >= min_readings[layer])
        outdata = pd.concat((outdata, all_targets.loc[idx, ['nztmx', 'nztmy', 'layer', 'h2o_elv_mean',
                                                            'weight', 'row', 'col', 'total_error_m']]))
    outdata.loc['BW23/0165'] = all_targets.loc['BW23/0165', ['nztmx', 'nztmy', 'layer', 'h2o_elv_mean',
                                                            'weight', 'row', 'col', 'total_error_m']]
    no_flow = smt.get_no_flow()
    for i in outdata.index:
        layer, row, col = outdata.loc[i, ['layer', 'row', 'col']].astype(int)
        if no_flow[layer, row, col] == 0:  # get rid of non-active wells
            outdata.loc[i, 'layer'] = np.nan
    outdata = outdata.dropna(subset=['layer', 'row', 'col'])
    # reviesd error for the NSMc as of 24/10/2017 take the total error and scale it between 1 and 10
    outdata.loc[:, 'error_2'] = (((outdata.loc[:, 'total_error_m'] - outdata.total_error_m.min()) /
                                  (outdata.total_error_m.max() - outdata.total_error_m.min())) * 9 + 1)
    outdata.loc[:, 'weight_2'] = 1 / outdata.loc[:, 'error_2']
    if return_missing:
        return outdata
    else:
        return outdata.loc[outdata.index != 'BW23/0165']


def get_seg_param_dict():
    seg_param_dict = {
        # small single point bits of the ashley
        8: 12,
        15: 18,
        22: 30,
        24: 30,
        28: 32,

        # I don't think these got implmented
        # link up small bits of the cust
        19: 21,
        37: 27,
        39: 27,
        41: 44,
        43: 44

    }
    return seg_param_dict


def save_dict_to_csv(path, indict, from_lab, to_lab):
    temp = pd.Series(indict, name=to_lab)
    temp.index.names = [from_lab]
    pd.DataFrame(temp).to_csv(path)


def check_non_head_targets():
    functions = [gen_constant_head_targets, gen_sfr_flux_target_array, gen_sfr_full_we_flux_target_array,
                 gen_sfr_flow_target_array, gen_drn_target_array]
    target_values = get_target_group_values()
    for f in functions:
        target_array, target_dict = f()
        if f == gen_constant_head_targets:
            for l in range(smt.layers):
                fig, ax = smt.plt_matrix(target_array[l], no_flow_layer=l)
                for j in target_dict.keys():
                    if target_values[target_dict[j]] is None:
                        temp_lab = target_values[target_dict[j]]
                    elif isinstance(target_values[target_dict[j]], str):
                        temp_lab = target_values[target_values[target_dict[j]]] / 86400
                    else:
                        temp_lab = target_values[target_dict[j]] / 86400
                    print('idx: {}, name: {}, value: {}'.format(j, target_dict[j], temp_lab))
                # plt.show(fig)
                plt.close(fig)
            continue
        for i in target_dict.keys():
            if target_values[target_dict[i]] is None:
                temp_lab = target_values[target_dict[i]]
            elif isinstance(target_values[target_dict[i]], str):
                temp_lab = target_values[target_values[target_dict[i]]] / 86400
            else:
                temp_lab = target_values[target_dict[i]] / 86400

            fig, ax = smt.plt_matrix(target_array,
                                     title='idx: {}, name: {}, value: {}'.format(i, target_dict[i], temp_lab),
                                     vmin=i - 0.5, vmax=i + 0.5)
            plt.show(fig)


if __name__ == '__main__':
    test = get_head_targets()

    import matplotlib.pyplot as plt
    temp = gen_sfr_flux_target_array()
    outdata = get_head_targets()
    outdata2 = get_head_targets(True)
    vert = get_vertical_gradient_targets()
    outdata.loc[:, 'test_error_1'] = (np.log10(outdata.loc[:, 'total_error_m']) + 1)
    outdata.loc[:, 'test_error_2'] = (((outdata.loc[:, 'total_error_m'] - outdata.total_error_m.min()) /
                                       (outdata.total_error_m.max() - outdata.total_error_m.min())) * 10 + 1)
    outdata.loc[:, 'test_error_3'] = outdata.loc[:, 'total_error_m']
    outdata.loc[outdata.test_error_3 < 1, 'test_error_3'] = 1
    outdata.loc[outdata.test_error_3 > 10, 'test_error_3'] = 10
    for var in ['total_error_m', 'weight', 'test_error_1', 'test_error_2', 'test_error_3','weight_2']:
        fig, ax = plt.subplots()
        data = []
        for k in range(smt.layers):
            data.append(outdata.loc[outdata.layer == k, var].values)
        ax.boxplot(data)
        ax.set_title(var)
    plt.show()
    print('done')
