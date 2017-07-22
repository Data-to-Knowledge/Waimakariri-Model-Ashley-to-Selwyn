# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 20/07/2017 12:55 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.m_packages.drn_packages import _get_drn_spd
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt, _get_constant_heads
import geopandas as gpd


# todo I need to get targets for everything from our data. most of the old targets are in the shp file.
# todo it might be worth popping a target dict in this script

def gen_drn_target_array():
    drn_data = _get_drn_spd(smt.reach_v, smt.wel_version)
    out_dict = {}
    out_array = np.zeros((smt.rows, smt.cols))
    for i, group in enumerate(set(drn_data['target_group']), 1):
        temp = drn_data.loc[drn_data.target_group == group]
        temp_array = smt.df_to_array(temp, 'k')
        out_array[np.isfinite(temp_array)] = i
        out_dict[i] = group

    return out_array, out_dict


def gen_sfr_flow_target_array():  # todo talk to cath/brioch about which to include, but can be done from existing shape files
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
                   14: 'sfo_e_wolf'}
    return target_array, num_to_name


def gen_sfr_full_we_flux_target_array():
    shp_path = '{}/m_ex_bd_inputs/shp/full_w_e_flux_targets.shp'.format(smt.sdp)
    target_array = smt.shape_file_to_model_array(shp_path, 'GRID_CODE', True)
    target_array[np.isnan(target_array)] = 0
    num_to_name = {1: 'sfx_w_all',
                   2: 'sfx_e_all'}

    return target_array, num_to_name

def gen_sfr_flux_target_array():  # todo talk to cath/brioch about which to include but can be done from existing shape files
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
                   115: 'sfx_1drn',
                   116: 'sfx_3drn',
                   117: 'sfx_4drn',
                   118: 'sfx_5drn',
                   119: 'sfx_2drn',
                   120: 'sfx_7drn'}

    return target_array, num_to_name


def gen_constant_head_targets():  # watch if we have constant heads in the sw boundary also note that this is 3d
    chbs = _get_constant_heads()
    shp_path = "{}/m_ex_bd_inputs/shp/coastal_head_target_zones.shp".format(smt.sdp)
    zones = np.repeat(smt.shape_file_to_model_array(shp_path, 'Id', alltouched=True)[np.newaxis, :, :], smt.layers,
                      axis=0)
    zones[np.isnan(chbs) | (chbs < -900)] = 0
    zones[np.isnan(zones)] = 0
    zone_data = gpd.read_file(shp_path)
    zone_data = zone_data.set_index('Id')
    zone_data = zone_data.loc[:, 'name'].to_dict()

    return zones, zone_data


def get_target_group_values():
    # Note that below are in m3/s and then converted to m3/day
    target_group_val = {'chb_ash': -5.25,  # 3.0 to 7.5 #todo talk over these
                        'chb_chch': -0.9,  # 0.3 to 1.5
                        'chb_cust': -0.3,  # 0.1 to 0.5
                        'chb_sely': 'sel_off',

                        # drains
                        'd_ash_car': None,
                        'd_ash_est': -1.2,  # 0.2 to 2
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
                        'sfo_c_pat': 0.0,
                        'sfo_c_skew': 1.75,
                        'sfo_c_swan': 0.85,
                        'sfo_c_tip': 0.0,
                        'sfo_c_tlks': 1.75,
                        'sfo_e_wolf': 0.0,

                        # surface water flux #from pervious shapefiles of targets
                        'sfx_a1_con': 5.20,
                        'sfx_a2_gol': 5.20,
                        'sfx_a3_tul': 5.20,
                        'sfx_a4_sh1': -0.40,
                        'sfx_c1_swa': -0.43,
                        'sfx_c2_mil': -0.42,
                        'sfx_e1_stf': 1.6,
                        'sfx_w1_cou': 2.8,
                        'sfx_w2_tom': 1.1,
                        'sfx_w3_ros': 2.2,
                        'sfx_w4_mcl': 5.7,
                        'sfx_w5_wat': -0.1,
                        'sfx_w6_sh1': -0.5,
                        'sfx_1drn': -0.23,
                        'sfx_2drn': -0.19,
                        'sfx_3drn': -0.18,
                        'sfx_4drn': -0.05,
                        'sfx_5drn': -0.05,
                        'sfx_7drn': -0.21,
                        'sfx_e_all': -1.99,
                        'sfx_w_all': -999999, #todo ask brioch what this target was

                        # groups
                        'sel_off': -11,  # 1.2-17
                        'chch_str': -10,  # 7.5 to 12.5 #todo it was -6 in the previous run think about this!!!
                        'sel_str': -9.8  # no range
                        }
    for key in target_group_val.keys():
        if isinstance(target_group_val[key], str) or target_group_val[key] is None:
            continue
        target_group_val[key] *= 86400 #convert numbers to m3/day

    return target_group_val


def get_vertical_gradient_targets():


    # load in vert targets
    vert_targets = pd.read_excel(env.sci(
        "Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/Vertical gradient targets updated_use.xlsx"),
                                 sheetname='data_for_python', index_col=0)

    # load in the row, col options from the bulk head targets sheet
    all_wells = pd.read_csv('{}/all_wells_row_col_layer.csv'.format(smt.sdp),index_col=0)
    vert_targets = pd.merge(vert_targets, all_wells,how='left',left_index=True,right_index=True)

    vert_targets.loc['M35/11937','GWL_RL'] = vert_targets.loc[['M35/11937','M35/10909'],'GWL_RL'].mean()
    vert_targets = vert_targets.drop(['M35/10909']) # this and above are in the same layer
    vert_targets.loc[:, 'weightings'] = None #todo weigthing

    outdata = vert_targets.loc[:,['NZTM_x','NZTM_y','layer','GWL_RL','weightings','row','col']].rename(columns={'NZTM_x':'x','NZTM_y':'y','GWL_RL':'obs','row':'i','col':'j'})

    #return a dataframe: lat, lon, layer, obs, weight?, i,j

    return outdata



if __name__ == '__main__':
    get_vertical_gradient_targets()
    zones, zone_data = gen_constant_head_targets()
    drn_array, drn_dict = gen_drn_target_array()
    flow_array, flow_dict = gen_sfr_flow_target_array()
    flux_array, flux_dict = gen_sfr_flux_target_array()
    print 'done'
