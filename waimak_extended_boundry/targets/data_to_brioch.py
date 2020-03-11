# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 1/09/2017 4:25 PM
"""

from __future__ import division
from waimak_extended_boundry.supporting_data_analysis.NSMC_inputs.data_to_cath import data_to_cath
from waimak_extended_boundry.m_packages.drn_packages import _get_drn_spd
from waimak_extended_boundry.m_packages.wel_packages import get_wel_spd
from waimak_extended_boundry.targets.gen_target_arrays import *
import os

def generate_all_data_for_brioch(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    head_targets = get_head_targets()
    head_targets.to_csv('{}/head_targets.csv'.format(outdir))

    vert_targets = get_vertical_gradient_targets()
    vert_targets.to_csv('{}/vertical_gradient_targets.csv'.format(outdir))

    target_group_val = get_target_group_values()
    save_dict_to_csv('{}/target_values.csv'.format(outdir), target_group_val, 'group_name', 'value')

    constant_head_targets, constant_head_zones = gen_constant_head_targets()
    for i in range(smt.layers):
        np.savetxt('{}/constant_head_target_0idlayer_{}.txt'.format(outdir, i), constant_head_targets[i])
    save_dict_to_csv('{}/constant_heads_dict.csv'.format(outdir), constant_head_zones, 'num_id', 'chb_id')

    futns = [gen_sfr_flux_target_array, gen_sfr_full_we_flux_target_array, gen_sfr_flow_target_array,
             gen_drn_target_array]
    names = ['sfr_flux', 'sfr_full_str_flux', 'sfr_flow', 'drn_flux']
    for f, n in zip(futns, names):
        temp_array, temp_dict = f()
        np.savetxt('{}/{}_array.txt'.format(outdir, n), temp_array)
        save_dict_to_csv('{}/{}_dict.csv'.format(outdir, n), temp_dict, 'num_id', '{}_id'.format(n))

    seg_dict = get_seg_param_dict()
    save_dict_to_csv('{}/segment_k_param_group.csv'.format(outdir), seg_dict, 'seg', 'use_k_of_seg')
    # well_ data
    well_data = get_wel_spd(smt.wel_version, True)
    well_data = well_data.loc[:, ['layer', 'row', 'col', 'flux', 'type']]
    well_data.to_csv('{}/well_data.csv'.format(outdir))

    drn_data = _get_drn_spd(smt.reach_v, smt.wel_version, True)
    drn_data.to_csv('{}/drain_data.csv'.format(outdir))
    waimak_springfed_targets = ['d_cam_mrsh',
                        'd_cam_revl',
                        'd_cam_yng',
                        'd_cour_nrd',
                        'd_emd_gard',
                        'd_kairaki',
                        'd_kuku_leg',
                        'd_nbk_mrsh',
                        'd_oho_btch',
                        'd_oho_jefs',
                        'd_oho_kpoi',
                        'd_oho_misc',
                        'd_oho_mlbk',
                        'd_oho_whit',
                        'd_salt_fct',
                        'd_salt_top',
                        'd_sbk_mrsh',
                        'd_sil_harp',
                        'd_sil_heyw',
                        'd_sil_ilnd',
                        'd_smiths',
                        'd_tar_gre',
                        'd_tar_stok']
    drn_data = smt.add_mxmy_to_df(drn_data)
    drn_data.loc[:,'has_target'] = np.in1d(drn_data.target_group, waimak_springfed_targets)
    drn_data.to_csv('{}/drain_data_for_shapefile.csv'.format(outdir))
    data_to_cath('{}'.format(outdir))
    print('done')


if __name__ == '__main__':
    temp = gen_constant_head_targets()
    temp = gen_sfr_flow_target_array()
    generate_all_data_for_brioch(r"C:\Users\MattH\Desktop\data_to_brioch_2017_10_24")