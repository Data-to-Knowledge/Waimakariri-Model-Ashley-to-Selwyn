# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 27/10/2017 10:25 AM
"""

from __future__ import division
from core import env
import netCDF4 as nc
from rrfextract import extractrrf, extractphisummary, extractoptphi
import numpy as np
import pandas as pd
import sys
import datetime
from read_pst import rd_pst_parameters, param_from_rec, extract_obs_opt_rei, extract_opt_priors
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import os

# rap up the NSMC parameters adn observations into a netcdf file
nsmc_dim = 7680 + 2
rch_dim = 46
layer_dim = 11
sfr_dim = 47
khv_dim = 178
well_dim = 1293
drn_dim = 37


def _add_simple_params(param, pst_param, prior_sd_data, postopt_sd_data, nc_file):
    simple_parameters = {'pump_c': {'units': 'none',
                                    'long_name': 'christchurch west melton pumping multiplier',
                                    'sd_type': 'lin',
                                    'opt_p': 1},
                         'pump_s': {'units': 'none',
                                    'long_name': 'selwyn pumping multiplier',
                                    'sd_type': 'lin',
                                    'opt_p': 1},
                         'pump_w': {'units': 'none',
                                    'long_name': 'waimakariri pumping multiplier',
                                    'sd_type': 'lin',
                                    'opt_p': 1},
                         'sriv': {'units': 'none',
                                  'long_name': 'selwyn river influx multiplier',
                                  'sd_type': 'lin',
                                  'opt_p': 1},
                         'n_race': {'units': 'none',
                                    'long_name': 'waimakariri race multiplier',
                                    'sd_type': 'lin',
                                    'opt_p': 1},
                         's_race': {'units': 'none',
                                    'long_name': 'selwyn race multiplier',
                                    'sd_type': 'lin',
                                    'opt_p': 1},
                         'nbndf': {'units': 'none',
                                   'long_name': 'northern boundary flux multiplier',
                                   'sd_type': 'lin',
                                   'opt_p': 1},
                         'top_e_flo': {'units': 'm3/day',
                                       'long_name': 'top of the eyre flow',
                                       'sd_type': 'log',
                                       'opt_p': 171936},
                         'mid_c_flo': {'units': 'm3/day',
                                       'long_name': 'mid cust (biwash) flow',
                                       'sd_type': 'log',
                                       'opt_p': 3456},
                         'top_c_flo': {'units': 'm3/day',
                                       'long_name': 'top of the cust flow',
                                       'sd_type': 'log',
                                       'opt_p': 14688},
                         'ulrzf': {'units': 'm3/day',
                                   'long_name': 'inland southwestern boundary flux',
                                   'sd_type': 'log',
                                   'opt_p': 259201.1
                                   },
                         'llrzf': {'units': 'm3/day',
                                   'long_name': 'coastal southwestern boundary flux',
                                   'sd_type': 'log',
                                   'opt_p': 181440.8},
                         'fkh_mult': {'units': 'none',
                                      'long_name': 'fault kh multiplier',
                                      'sd_type': 'log',
                                      'opt_p': 1},
                         'fkv_mult': {'units': 'none',
                                      'long_name': 'fault kv multiplier',
                                      'sd_type': 'log',
                                      'opt_p': 1}}

    for key in simple_parameters.keys():
        temp = nc_file.createVariable(key, 'f8', ('nsmc_num',), fill_value=np.nan, zlib=False)
        temp.setncatts({'units': simple_parameters[key]['units'],
                        'long_name': simple_parameters[key]['long_name'],
                        'missing_value': np.nan,
                        'initial': pst_param.loc[key, 'initial'],
                        'lower': pst_param.loc[key, 'lower'],
                        'upper': pst_param.loc[key, 'upper'],
                        'vtype': 'param',
                        'sd_type': simple_parameters[key]['sd_type'],
                        'p_sd': prior_sd_data.loc[key, 'sd'],
                        'j_sd': postopt_sd_data.loc[key, 'sd'] ** 0.5,
                        'opt_p': simple_parameters[key]['opt_p']
                        })
        adder = 0
        if key == 'llrzf':
            adder = -181441  # offset in pest
        elif key == 'ulrzf':
            adder = -1

        temp[:] = param.loc[key].values + adder


def _add_rch_params(param, rch_ppt_tpl, pst_param, prior_sd_data, postopt_sd_data, nc_file):
    # rch ppts
    rch_meta = pd.read_table(rch_ppt_tpl, skiprows=1, names=['x', 'y', 'group', 'ignore'], delim_whitespace=True)
    rch_ppt_ids = ['rch_ppt_{:02d}'.format(e) for e in range(46)]
    rch_pts = nc_file.createVariable('rch_ppt', str, ('rch_ppt',), zlib=False)
    rch_pts.setncatts({'units': 'none',
                       'long_name': 'recharge pilot point identifier',
                       'comments': 'this is a unique identifier',
                       'vtype': 'dim'})
    for i, val in enumerate(rch_ppt_ids):
        rch_pts[i] = val

    rch_x = nc_file.createVariable('rch_ppt_x', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_x.setncatts({'units': 'nztmx',
                     'long_name': 'recharge pilot point longitude',
                     'missing_value': np.nan,
                     'vtype': 'meta'})
    rch_x[:] = rch_meta.loc[:, 'x'].values

    rch_y = nc_file.createVariable('rch_ppt_y', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_y.setncatts({'units': 'nztmy',
                     'long_name': 'recharge pilot point latitude',
                     'missing_value': np.nan,
                     'vtype': 'meta'})
    rch_y[:] = rch_meta.loc[:, 'y'].values

    rch_opt_p = nc_file.createVariable('rch_ppt_opt_p', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_opt_p.setncatts({'units': 'nztmy',
                         'long_name': 'recharge pilot point optimisation prior',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    rch_opt_p[:] = np.ones((len(rch_ppt_ids)))

    # initials uppers lowers
    rch_lower = nc_file.createVariable('rch_ppt_lower', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_lower.setncatts({'units': 'none',
                         'long_name': 'recharge pilot point lower bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    rch_lower[:] = pst_param.loc[rch_ppt_ids, 'lower'].values

    rch_upper = nc_file.createVariable('rch_ppt_upper', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_upper.setncatts({'units': 'none',
                         'long_name': 'recharge pilot point upper bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    rch_upper[:] = pst_param.loc[rch_ppt_ids, 'upper'].values

    rch_initial = nc_file.createVariable('rch_ppt_initial', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_initial.setncatts({'units': 'none',
                           'long_name': 'recharge pilot point initial value',
                           'missing_value': np.nan,
                           'vtype': 'meta'})
    rch_initial[:] = pst_param.loc[rch_ppt_ids, 'initial'].values

    rch_sd_prior = nc_file.createVariable('rch_ppt_p_sd', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_sd_prior.setncatts({'units': 'none',
                            'long_name': 'recharge pilot point prior standard deviation',
                            'missing_value': np.nan,
                            'vtype': 'meta',
                            'sd_type': 'lin'})
    rch_sd_prior[:] = prior_sd_data.loc[rch_ppt_ids, 'sd'].values

    rch_sd_post = nc_file.createVariable('rch_ppt_j_sd', 'f8', ('rch_ppt',), fill_value=np.nan, zlib=False)
    rch_sd_post.setncatts({'units': 'none',
                           'long_name': 'recharge pilot point post sensitivity standard deviation',
                           'missing_value': np.nan,
                           'vtype': 'meta',
                           'sd_type': 'lin'})
    rch_sd_post[:] = postopt_sd_data.loc[rch_ppt_ids, 'sd'].values ** 0.5

    rch_group = nc_file.createVariable('rch_ppt_group', 'i4', ('rch_ppt',), fill_value=-9, zlib=False)
    rch_group.setncatts({'flag_values': [1, 2, 3, 4],
                         'flag_meanings': 'dryland confined selwyn_irr waimak_irr',
                         'long_name': 'recharge pilot point groups',
                         'missing_value': -9,
                         'vtype': 'meta'})
    rch_group[:] = rch_meta.loc[:, 'group'].values

    rch_mult = nc_file.createVariable('rch_mult', 'f8', ('nsmc_num', 'rch_ppt'), fill_value=np.nan, zlib=False)
    rch_mult.setncatts({'units': 'none',
                        'long_name': 'recharge multipliers',
                        'missing_value': np.nan,
                        'vtype': 'param'})
    temp_data = np.zeros((nsmc_dim, rch_dim)) * np.nan
    for i, key in enumerate(rch_ppt_ids):
        temp_data[:, i] = param.loc[key].values

    rch_mult[:] = temp_data


def _add_sfr_cond(param, pst_param, prior_sd_data, postopt_sd_data, nc_file):
    # stream hconds
    hcond_sites = ['hcond1', 'hcond2', 'hcond3', 'hcond4', 'hcond5', 'hcond6', 'hcond7', 'hcond8', 'hcond9', 'hcond10',
                   'hcond11', 'hcond12', 'hcond13', 'hcond14', 'hcond15', 'hcond16', 'hcond17', 'hcond18', 'hcond19',
                   'hcond20', 'hcond21', 'hcond22', 'hcond23', 'hcond24', 'hcond25', 'hcond26', 'hcond27', 'hcond28',
                   'hcond29', 'hcond30', 'hcond31', 'hcond32', 'hcond33', 'hcond34', 'hcond34x', 'hcond35', 'hcond35x',
                   'hcond36', 'hcond37', 'hcond38', 'hcond39', 'hcond40', 'hcond41', 'hcond42', 'hcond43', 'hcond44',
                   'hcond44x']

    sfr_cond = nc_file.createVariable('sfr_cond', str, ('sfr_cond',), zlib=False)
    sfr_cond.setncatts({'units': 'none',
                        'long_name': 'sfr conductance identifier',
                        'vtype': 'dim'})
    for i, val in enumerate(hcond_sites):
        sfr_cond[i] = val

    # I could add river that the segment relates to if I have time (I wont)

    # add upper and lower bounds
    sfr_lower = nc_file.createVariable('sfr_lower', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_lower.setncatts({'units': 'none',
                         'long_name': 'sfr cond lower bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    sfr_lower[:] = pst_param.loc[hcond_sites, 'lower'].values

    sfr_upper = nc_file.createVariable('sfr_upper', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_upper.setncatts({'units': 'none',
                         'long_name': 'sfr cond upper bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    sfr_upper[:] = pst_param.loc[hcond_sites, 'upper'].values

    sfr_initial = nc_file.createVariable('sfr_initial', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_initial.setncatts({'units': 'none',
                           'long_name': 'sfr cond initial value',
                           'missing_value': np.nan,
                           'vtype': 'meta'})
    sfr_initial[:] = pst_param.loc[hcond_sites, 'initial'].values

    sfr_prior_sd = nc_file.createVariable('sfr_p_sd', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_prior_sd.setncatts({'units': 'none',
                            'long_name': 'sfr cond prior standard deviation',
                            'missing_value': np.nan,
                            'vtype': 'meta',
                            'sd_type': 'log'
                            })
    sfr_prior_sd[:] = prior_sd_data.loc[hcond_sites, 'sd'].values

    sfr_opt_p = nc_file.createVariable('sfr_opt_p', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_opt_p.setncatts({'units': 'none',
                         'long_name': 'sfr cond optimisation prior ',
                         'missing_value': np.nan,
                         'vtype': 'meta',
                         'sd_type': 'log'
                         })
    sfr_opt_p[:] = np.ones(len(hcond_sites)).fill(10)

    sfr_post_sd = nc_file.createVariable('sfr_j_sd', 'f8', ('sfr_cond',), fill_value=np.nan, zlib=False)
    sfr_post_sd.setncatts({'units': 'none',
                           'long_name': 'sfr cond post sensitivity matrix standard deviation',
                           'missing_value': np.nan,
                           'vtype': 'meta',
                           'sd_type': 'log'
                           })
    sfr_post_sd[:] = postopt_sd_data.loc[hcond_sites, 'sd'].values ** 0.5

    sfr_cond_val = nc_file.createVariable('sfr_cond_val', 'f8', ('nsmc_num', 'sfr_cond'), fill_value=np.nan, zlib=False)
    sfr_cond_val.setncatts({'units': 'm/day',
                            'long_name': 'sfr conductance at points',
                            'missing_value': np.nan,
                            'vtype': 'param'})

    temp_data = np.zeros((nsmc_dim, sfr_dim)) * np.nan
    for i, key in enumerate(hcond_sites):
        temp_data[:, i] = param.loc[key].values
    sfr_cond_val[:] = temp_data


def _add_drain_cond(param, pst_param, prior_sd_data, postopt_sd_data, nc_file):
    drns = ['d_salt_fct', 'd_salt_top', 'd_kuku_leg', 'd_ash_s', 'd_tar_stok', 'd_cam_mrsh', 'd_tar_gre', 'd_nbk_mrsh',
            'd_sbk_mrsh', 'd_cam_yng', 'd_cam_revl', 'd_kairaki', 'd_oho_mlbk', 'd_oho_whit', 'd_oho_jefs',
            'd_oho_kpoi', 'd_sil_ilnd', 'd_oho_misc', 'd_oho_btch', 'd_cour_nrd', 'd_dwaimak', 'd_emd_gard',
            'd_bul_styx', 'd_sil_heyw', 'd_smiths', 'd_sil_harp', 'd_bul_avon', 'd_ash_c', 'd_cust_c', 'd_chch_c',
            'd_waihora', 'd_ulin_c', 'd_dlin_c', 'd_usel_c', 'd_dsel_c', 'd_uwaimak', 'd_ash_est']
    drn_pts = nc_file.createVariable('drns', str, ('drns',), zlib=False)
    drn_pts.setncatts({'units': 'none',
                       'long_name': 'drain_segment names',
                       'vtype': 'dim'})
    for i, val in enumerate(drns):
        drn_pts[i] = val

    drn_lower = nc_file.createVariable('drn_lower', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_lower.setncatts({'units': 'none',
                         'long_name': 'drain cond lower bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    drn_lower[:] = pst_param.loc[drns, 'lower'].values

    drn_upper = nc_file.createVariable('drn_upper', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_upper.setncatts({'units': 'none',
                         'long_name': 'drain cond upper bound',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    drn_upper[:] = pst_param.loc[drns, 'upper'].values

    drn_initial = nc_file.createVariable('drn_initial', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_initial.setncatts({'units': 'none',
                           'long_name': 'drain cond initial value',
                           'missing_value': np.nan,
                           'vtype': 'meta'})
    drn_initial[:] = pst_param.loc[drns, 'initial'].values

    drn_prior_sd = nc_file.createVariable('drn_p_sd', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_prior_sd.setncatts({'units': 'none',
                            'long_name': 'drain cond prior standard deviation',
                            'missing_value': np.nan,
                            'vtype': 'meta',
                            'sd_type': 'log'})
    drn_prior_sd[:] = prior_sd_data.loc[drns, 'sd'].values

    drn_opt_p = nc_file.createVariable('drn_opt_p', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_opt_p.setncatts({'units': 'none',
                         'long_name': 'drain cond optimisation prior',
                         'missing_value': np.nan,
                         'vtype': 'meta',
                         'sd_type': 'log'})
    drn_opt_p[:] = np.ones(len(drns)).fill(1680)

    drn_post_sd = nc_file.createVariable('drn_j_sd', 'f8', ('drns',), fill_value=np.nan, zlib=False)
    drn_post_sd.setncatts({'units': 'none',
                           'long_name': 'drain cond post sensitivity matrix standard deviation',
                           'missing_value': np.nan,
                           'vtype': 'meta',
                           'sd_type': 'log'})
    drn_post_sd[:] = postopt_sd_data.loc[drns, 'sd'].values ** 0.5

    drn_cond = nc_file.createVariable('drn_cond', 'f8', ('nsmc_num', 'drns'), fill_value=np.nan, zlib=False)
    drn_cond.setncatts({'units': 'm3/day',
                        'long_name': 'drain conductance',
                        'missing_value': np.nan,
                        'vtype': 'param'})
    drn_cond[:] = param.loc[drns].transpose().values


def _add_kv_kh(param, kh_kv_ppt_file, pst_param, prior_ksds_dir, postopt_sd_data, opt_pst, nc_file):
    opt_priors = extract_opt_priors(opt_pst)

    # kh kv pilot points dimensions (nsmc_num, ppt, layer)
    ppt_ids = ['pp005602', 'pp005630', 'pp005798', 'pp010726', 'pp010754', 'pp010810', 'pp010894', 'pp010922',
               'pp015738', 'pp015794', 'pp015822', 'pp015850', 'pp015878', 'pp015906', 'pp015934', 'pp015962',
               'pp015990', 'pp016018', 'pp020834', 'pp020862', 'pp020890', 'pp020918', 'pp020946', 'pp020974',
               'pp021002', 'pp021030', 'pp021058', 'pp021086', 'pp021114', 'pp021142', 'pp025930', 'pp025958',
               'pp025986', 'pp026014', 'pp026042', 'pp026070', 'pp026098', 'pp026126', 'pp026154', 'pp026182',
               'pp026210', 'pp026238', 'pp031054', 'pp031082', 'pp031110', 'pp031138', 'pp031166', 'pp031194',
               'pp031222', 'pp031250', 'pp031278', 'pp031306', 'pp031334', 'pp031362', 'pp036150', 'pp036178',
               'pp036206', 'pp036234', 'pp036262', 'pp036290', 'pp036318', 'pp036346', 'pp036374', 'pp036402',
               'pp036430', 'pp036458', 'pp041302', 'pp041330', 'pp041358', 'pp041386', 'pp041414', 'pp041442',
               'pp041470', 'pp041498', 'pp041526', 'pp041554', 'pp041582', 'pp046398', 'pp046426', 'pp046454',
               'pp046482', 'pp046510', 'pp046538', 'pp046566', 'pp046594', 'pp046622', 'pp046650', 'pp046678',
               'pp051522', 'pp051550', 'pp051578', 'pp051606', 'pp051634', 'pp051662', 'pp051690', 'pp051718',
               'pp051746', 'pp051774', 'pp051802', 'pp056618', 'pp056646', 'pp056674', 'pp056702', 'pp056730',
               'pp056758', 'pp056786', 'pp056814', 'pp056842', 'pp056870', 'pp056898', 'pp061714', 'pp061742',
               'pp061770', 'pp061798', 'pp061826', 'pp061854', 'pp061882', 'pp061910', 'pp061938', 'pp061966',
               'pp061994', 'pp062022', 'pp066838', 'pp066866', 'pp066894', 'pp066922', 'pp066950', 'pp066978',
               'pp067006', 'pp067034', 'pp067062', 'pp067090', 'pp067118', 'pp071962', 'pp071990', 'pp072018',
               'pp072046', 'pp072074', 'pp072102', 'pp072130', 'pp072158', 'pp072186', 'pp079628', 'pp079668',
               'pp079708', 'pp079748', 'pp079788', 'pp079828', 'pp086948', 'pp086988', 'pp087028', 'pp087068',
               'pp087108', 'pp094268', 'pp094308', 'pp094348', 'pp094388', 'pp094428', 'pp101588', 'pp101628',
               'pp101668', 'pp101708', 'pp108908', 'pp108948', 'pp108988', 'pp109028', 'pp116228', 'pp116268',
               'pp116308', 'pp116348', 'pp123548', 'pp123588', 'pp123628', 'pp123668', 'pp130868', 'pp130908',
               'pp130948', 'pp130988']  # from layer 1 I assume this is all of them
    ppt_meta = pd.read_table(kh_kv_ppt_file, skiprows=1, names=['x', 'y', 'group', 'ignore'], delim_whitespace=True)

    ppts = nc_file.createVariable('khv_ppt', str, ('khv_ppt',), zlib=False)
    ppts.setncatts({'units': 'none',
                    'long_name': 'kv and kh pilot point id',
                    'vtype': 'dim'})
    for i, val in enumerate(ppt_ids):
        ppts[i] = val

    # X
    pptsx = nc_file.createVariable('khv_ppt_x', 'f8', ('khv_ppt',), fill_value=np.nan, zlib=False)
    pptsx.setncatts({'units': 'nztm',
                     'long_name': 'kv and kh pilot point longitude',
                     'missing_value': np.nan,
                     'vtype': 'meta'})
    pptsx[:] = ppt_meta.loc[ppt_ids, 'x'].values

    # Y
    pptsy = nc_file.createVariable('khv_ppt_y', 'f8', ('khv_ppt',), fill_value=np.nan, zlib=False)
    pptsy.setncatts({'units': 'nztm',
                     'long_name': 'kv and kh pilot point latitude',
                     'missing_value': np.nan,
                     'vtype': 'meta'})
    pptsy[:] = ppt_meta.loc[ppt_ids, 'y'].values

    # kv upper lower initial
    temp_upper = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_lower = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_initial = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_prior_sd = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_post_sd = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_opt_p = np.zeros((layer_dim, khv_dim)) * np.nan
    for i in range(0, layer_dim):
        prior_ksd = pd.read_table(os.path.join(prior_ksds_dir, '{:02}_cov_v_07.txt'.format(i + 1)),
                                  skiprows=1, index_col=0, names=['sd'], delim_whitespace=True)
        for k, key in enumerate(ppt_ids):
            try:
                temp_upper[i, k] = pst_param.loc['{}_v{}'.format(key, i + 1), 'upper']
                temp_lower[i, k] = pst_param.loc['{}_v{}'.format(key, i + 1), 'lower']
                temp_initial[i, k] = pst_param.loc['{}_v{}'.format(key, i + 1), 'initial']
                temp_prior_sd[i, k] = prior_ksd.loc['{}_v{}'.format(key, i + 1), 'sd']
                temp_post_sd[i, k] = postopt_sd_data.loc['{}_v{}'.format(key, i + 1), 'sd']
                temp_opt_p[i, k] = opt_priors.loc['{}_v{}'.format(key, i + 1)]
            except KeyError:
                pass

    pptlower = nc_file.createVariable('kv_lower', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptlower.setncatts({'units': 'm/day',
                        'long_name': 'kv pilot point lower bound',
                        'missing_value': np.nan,
                        'vtype': 'meta'})
    pptlower[:] = temp_lower

    pptupper = nc_file.createVariable('kv_upper', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptupper.setncatts({'units': 'm/day',
                        'long_name': 'kv pilot point upper bound',
                        'missing_value': np.nan,
                        'vtype': 'meta'})
    pptupper[:] = temp_upper

    pptinitial = nc_file.createVariable('kv_initial', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptinitial.setncatts({'units': 'm/day',
                          'long_name': 'kv pilot point initial value',
                          'missing_value': np.nan,
                          'vtype': 'meta'})
    pptinitial[:] = temp_initial

    pptprior_sd = nc_file.createVariable('kv_p_sd', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptprior_sd.setncatts({'units': 'm/day',
                           'long_name': 'kv pilot point prior standard deviation',
                           'missing_value': np.nan,
                           'vtype': 'meta',
                           'sd_type': 'log'})
    pptprior_sd[:] = temp_prior_sd ** 0.5

    pptoptp = nc_file.createVariable('kv_opt_p', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptoptp.setncatts({'units': 'm/day',
                       'long_name': 'kv pilot point optimisation prior',
                       'missing_value': np.nan,
                       'vtype': 'meta',
                       'sd_type': 'log'})

    pptoptp[:] = temp_opt_p

    pptpost_sd = nc_file.createVariable('kv_j_sd', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptpost_sd.setncatts({'units': 'm/day',
                          'long_name': 'kv pilot point post sensitivity matrix standard deviation',
                          'missing_value': np.nan,
                          'vtype': 'meta',
                          'sd_type': 'log'})
    pptpost_sd[:] = temp_post_sd ** 0.5

    # kh upper lower initial
    temp_upper = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_lower = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_initial = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_prior_sd = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_post_sd = np.zeros((layer_dim, khv_dim)) * np.nan
    temp_opt_p = np.zeros((layer_dim, khv_dim)) * np.nan
    for i in range(0, layer_dim):
        prior_ksd = pd.read_table(os.path.join(prior_ksds_dir, '{:02}_cov_05.txt'.format(i + 1)),
                                  skiprows=1, index_col=0, names=['sd'], delim_whitespace=True)

        for k, key in enumerate(ppt_ids):
            try:
                temp_upper[i, k] = pst_param.loc['{}_h{}'.format(key, i + 1), 'upper']
                temp_lower[i, k] = pst_param.loc['{}_h{}'.format(key, i + 1), 'lower']
                temp_initial[i, k] = pst_param.loc['{}_h{}'.format(key, i + 1), 'initial']
                temp_prior_sd[i, k] = prior_ksd.loc['{}_h{}'.format(key, i + 1), 'sd']
                temp_post_sd[i, k] = postopt_sd_data.loc['{}_h{}'.format(key, i + 1), 'sd']
                temp_opt_p[i, k] = opt_priors.loc['{}_h{}'.format(key, i + 1)]

            except KeyError:
                pass

    pptlower = nc_file.createVariable('kh_lower', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptlower.setncatts({'units': 'm/day',
                        'long_name': 'kh pilot point lower bound',
                        'missing_value': np.nan,
                        'vtype': 'meta'})
    pptlower[:] = temp_lower

    pptupper = nc_file.createVariable('kh_upper', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptupper.setncatts({'units': 'm/day',
                        'long_name': 'kh pilot point upper bound',
                        'missing_value': np.nan,
                        'vtype': 'meta'})
    pptupper[:] = temp_upper

    pptinitial = nc_file.createVariable('kh_initial', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptinitial.setncatts({'units': 'm/day',
                          'long_name': 'kh pilot point initial value',
                          'missing_value': np.nan,
                          'vtype': 'meta'})
    pptinitial[:] = temp_initial

    pptoptp = nc_file.createVariable('kh_opt_p', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptoptp.setncatts({'units': 'm/day',
                       'long_name': 'kh pilot point optimisation prior',
                       'missing_value': np.nan,
                       'vtype': 'meta',
                       'sd_type': 'log'})

    pptoptp[:] = temp_opt_p

    pptprior_sd = nc_file.createVariable('kh_p_sd', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptprior_sd.setncatts({'units': 'm/day',
                           'long_name': 'kh pilot point prior standard deviation',
                           'missing_value': np.nan,
                           'vtype': 'meta',
                           'sd_type': 'log'})
    pptprior_sd[:] = temp_prior_sd ** 0.5

    pptpost_sd = nc_file.createVariable('kh_j_sd', 'f8', ('layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    pptpost_sd.setncatts({'units': 'm/day',
                          'long_name': 'kh pilot point post sensitivity matrix standard deviation',
                          'missing_value': np.nan,
                          'vtype': 'meta',
                          'sd_type': 'log'})
    pptpost_sd[:] = temp_post_sd ** 0.5

    # kv
    kv = nc_file.createVariable('kv', 'f8', ('nsmc_num', 'layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    kv.setncatts({'units': 'm/day',
                  'long_name': 'vertical conductivity',
                  'missing_value': np.nan,
                  'vtype': 'param'})
    temp_data = np.zeros((nsmc_dim, layer_dim, khv_dim)) * np.nan
    for i in range(0, layer_dim):
        for k, key in enumerate(ppt_ids):
            try:
                temp_data[:, i, k] = param.loc['{}_v{}'.format(key, i + 1)]
            except KeyError:
                pass
    kv[:] = temp_data

    # kh
    kh = nc_file.createVariable('kh', 'f8', ('nsmc_num', 'layer', 'khv_ppt'), fill_value=np.nan, zlib=False)
    kh.setncatts({'units': 'm/day',
                  'long_name': 'horizontal conductivity',
                  'missing_value': np.nan,
                  'vtype': 'param'})
    temp_data = np.zeros((nsmc_dim, layer_dim, khv_dim)) * np.nan
    for i in range(0, layer_dim):
        for k, key in enumerate(ppt_ids):
            try:
                temp_data[:, i, k] = param.loc['{}_h{}'.format(key, i + 1)]
            except KeyError:
                pass
    kh[:] = temp_data


def _add_well_obs(obs_file, rei_file, nc_file):
    other_obs = ['brnthl_4_1', 'eyrftn_2_1', 'eyrftn_6_2', 'eyrftn_6_1', 'wdnd_4_2', 'wdnd_8_4', 'wdnd_8_2', 'peg_8_7',
                 'peg_9_8', 'peg_10_9', 'peg_10_7', 'eyrftm_2_1', 'eyrftl_7_2', 'oxfds_6_2', 'oxfdnr_3_1', 'oxfdnr_4_3',
                 'oxfdnr_4_1', 'chch_4_2', 'chb_ash', 'chb_chch', 'chb_chchz', 'chb_cust', 'chb_sely', 'sel_off',
                 'chch_str', 'd_ash_c', 'd_ash_est', 'd_ash_s', 'd_bul_avon', 'd_bul_styx', 'd_cam_mrsh', 'd_cam_revl',
                 'd_cam_revlz', 'd_cam_s', 'd_cam_yng', 'd_chch_c', 'd_cour_nrd', 'd_court_s', 'd_cust_c', 'sel_str',
                 'd_dlin_c', 'd_dsel_c', 'd_dwaimak', 'd_ect', 'd_emd_gard', 'd_kaiapo_s', 'd_kairaki', 'd_kuku_leg',
                 'd_nbk_mrsh', 'd_oho_btch', 'd_oho_jefs', 'd_oho_kpoi', 'd_oho_misc', 'd_oho_miscz', 'd_oho_mlbk',
                 'd_oho_whit', 'd_salt_fct', 'd_salt_s', 'd_salt_top', 'd_sbk_mrsh', 'd_sil_harp', 'd_sil_heyw',
                 'd_sil_ilnd', 'd_smiths', 'd_tar_gre', 'd_tar_stok', 'd_taran_s', 'd_ulin_c', 'd_usel_c', 'd_uwaimak',
                 'd_waihora', 'd_waikuk_s', 'sfx_2adrn', 'sfx_3drn', 'sfx_4drn', 'sfx_5drn', 'sfx_6drn', 'sfx_7drn',
                 'mid_ash_g', 'sfx_a1_con', 'sfx_a2_gol', 'sfx_a3_tul', 'sfx_a4_sh1', 'sfx_c1_swa', 'sfx_c2_mil',
                 'sfx_cd_all', 'sfx_e1_stf', 'sfx_e_all', 'sfx_w1_cou', 'sfx_w2_tom', 'sfx_w3_ros', 'sfx_w4_mcl',
                 'sfx_w5_wat', 'sfx_w6_sh1', 'sfx_w_all', 'sfo_1drn', 'sfo_2drn', 'sfo_3drn', 'sfo_4drn', 'sfo_5drn',
                 'sfo_7drn', 'sfo_c_benn', 'sfo_c_oxf', 'sfo_c_pat', 'sfo_c_skew', 'sfo_c_swan', 'sfo_c_tip',
                 'sfo_c_tlks', 'sfo_e_down', 'sfo_e_downz', 'sfo_e_poyn', 'sfo_e_poynz', 'sfo_e_seyr', 'sfo_e_seyrz',
                 'sfo_e_wolf', 'sfo_e_wolz']
    meta = pd.read_table(rei_file, skiprows=4, index_col=0, delim_whitespace=True)
    well_names = list(set(obs_file.index) - set(other_obs))
    well_names_cap = [e.upper().replace('_', '/') for e in well_names]

    # names
    well_nm = nc_file.createVariable('well_name', str, ('well_name',), zlib=False)
    well_nm.setncatts({'units': 'none',
                       'long_name': 'identifier for head targets',
                       'vtype': 'dim'})
    for i, val in enumerate(well_names_cap):
        well_nm[i] = val

    # obs
    obs = nc_file.createVariable('well_obs', 'f8', ('nsmc_num', 'well_name'), fill_value=np.nan, zlib=False)
    obs.setncatts({'units': 'm',
                   'long_name': 'model observation for head targets',
                   'missing_value': np.nan,
                   'vtype': 'obs'})
    obs[:] = obs_file.loc[well_names].transpose().values

    # target values
    tar = nc_file.createVariable('well_target', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    tar.setncatts({'units': 'm',
                   'long_name': 'head targets value',
                   'missing_value': np.nan,
                   'vtype': 'meta'})
    tar[:] = meta.loc[well_names, 'Measured'].values

    # weight
    weight = nc_file.createVariable('well_weight', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    weight.setncatts({'units': 'none',
                      'long_name': 'NSMC weight for head targets',
                      'missing_value': np.nan,
                      'vtype': 'meta'})
    weight[:] = meta.loc[well_names, 'Weight'].values

    all_wells = get_all_well_row_col()

    # x
    x = nc_file.createVariable('well_x', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    x.setncatts({'units': 'nztm',
                 'long_name': 'head targets longitude',
                 'missing_value': np.nan,
                 'vtype': 'meta'})
    x[:] = all_wells.loc[well_names_cap, 'nztmx'].values

    # y
    y = nc_file.createVariable('well_y', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    y.setncatts({'units': 'nztm',
                 'long_name': 'head targets latitude',
                 'missing_value': np.nan,
                 'vtype': 'meta'})
    y[:] = all_wells.loc[well_names_cap, 'nztmy'].values

    # depth
    depth = nc_file.createVariable('well_depth', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    depth.setncatts({'units': 'm',
                     'long_name': 'head targets depth',
                     'missing_value': np.nan,
                     'vtype': 'meta'})
    depth[:] = all_wells.loc[well_names_cap, 'depth'].values

    # midscreen_elv
    midscreen = nc_file.createVariable('well_midscreen', 'f8', ('well_name',), fill_value=np.nan, zlib=False)
    midscreen.setncatts({'units': 'm',
                         'long_name': 'head targets midscreen elevation',
                         'missing_value': np.nan,
                         'vtype': 'meta'})
    midscreen[:] = all_wells.loc[well_names_cap, 'mid_screen_elv'].values


def _add_other_obs(obs_file, rei_file, nc_file):
    other_obs = ['brnthl_4_1', 'eyrftn_2_1', 'eyrftn_6_2', 'eyrftn_6_1', 'wdnd_4_2', 'wdnd_8_4', 'wdnd_8_2', 'peg_8_7',
                 'peg_9_8', 'peg_10_9', 'peg_10_7', 'eyrftm_2_1', 'eyrftl_7_2', 'oxfds_6_2', 'oxfdnr_3_1', 'oxfdnr_4_3',
                 'oxfdnr_4_1', 'chch_4_2', 'chb_ash', 'chb_chch', 'chb_chchz', 'chb_cust', 'chb_sely', 'sel_off',
                 'chch_str', 'd_ash_c', 'd_ash_est', 'd_ash_s', 'd_bul_avon', 'd_bul_styx', 'd_cam_mrsh', 'd_cam_revl',
                 'd_cam_revlz', 'd_cam_s', 'd_cam_yng', 'd_chch_c', 'd_cour_nrd', 'd_court_s', 'd_cust_c', 'sel_str',
                 'd_dlin_c', 'd_dsel_c', 'd_dwaimak', 'd_ect', 'd_emd_gard', 'd_kaiapo_s', 'd_kairaki', 'd_kuku_leg',
                 'd_nbk_mrsh', 'd_oho_btch', 'd_oho_jefs', 'd_oho_kpoi', 'd_oho_misc', 'd_oho_miscz', 'd_oho_mlbk',
                 'd_oho_whit', 'd_salt_fct', 'd_salt_s', 'd_salt_top', 'd_sbk_mrsh', 'd_sil_harp', 'd_sil_heyw',
                 'd_sil_ilnd', 'd_smiths', 'd_tar_gre', 'd_tar_stok', 'd_taran_s', 'd_ulin_c', 'd_usel_c', 'd_uwaimak',
                 'd_waihora', 'd_waikuk_s', 'sfx_2adrn', 'sfx_3drn', 'sfx_4drn', 'sfx_5drn', 'sfx_6drn', 'sfx_7drn',
                 'mid_ash_g', 'sfx_a1_con', 'sfx_a2_gol', 'sfx_a3_tul', 'sfx_a4_sh1', 'sfx_c1_swa', 'sfx_c2_mil',
                 'sfx_cd_all', 'sfx_e1_stf', 'sfx_e_all', 'sfx_w1_cou', 'sfx_w2_tom', 'sfx_w3_ros', 'sfx_w4_mcl',
                 'sfx_w5_wat', 'sfx_w6_sh1', 'sfx_w_all', 'sfo_1drn', 'sfo_2drn', 'sfo_3drn', 'sfo_4drn', 'sfo_5drn',
                 'sfo_7drn', 'sfo_c_benn', 'sfo_c_oxf', 'sfo_c_pat', 'sfo_c_skew', 'sfo_c_swan', 'sfo_c_tip',
                 'sfo_c_tlks', 'sfo_e_down', 'sfo_e_downz', 'sfo_e_poyn', 'sfo_e_poynz', 'sfo_e_seyr', 'sfo_e_seyrz',
                 'sfo_e_wolf', 'sfo_e_wolz']

    verts = ['brnthl_4_1', 'eyrftn_2_1', 'eyrftn_6_2', 'eyrftn_6_1', 'wdnd_4_2', 'wdnd_8_4', 'wdnd_8_2', 'peg_8_7',
             'peg_9_8', 'peg_10_9', 'peg_10_7', 'eyrftm_2_1', 'eyrftl_7_2', 'oxfds_6_2', 'oxfdnr_3_1', 'oxfdnr_4_3',
             'oxfdnr_4_1', 'chch_4_2']

    meta = pd.read_table(rei_file, skiprows=4, index_col=0, delim_whitespace=True)

    for obs in other_obs:
        units = 'm3/day'
        if obs in verts:
            units = 'm'
            long_name = 'verticle gradient target {}'.format(obs)
        elif 'mid_ash_g' == obs:
            long_name = 'midashley stream-model flux group'
        elif 'sel_str' == obs:
            long_name = 'selwyn springfed stream flux'
        elif 'chch_str' == obs:
            long_name = 'chch springfed stream flux'
        elif 'sel_off' == obs:
            long_name = 'selwyn offshore flux'
        elif 'd_' in obs:
            long_name = 'drain flux {}'.format(obs)
        elif 'sfx_' in obs:
            long_name = 'sfr stream-model flux {}'.format(obs)
        elif 'chb_' in obs:
            long_name = 'constant head/offshore {}'.format(obs)
        elif 'sfo_' in obs:
            long_name = 'sfr in stream flow {}'.format(obs)
        else:
            raise ValueError('{} shouldnt get here'.format(obs))

        temp = nc_file.createVariable(obs, 'f8', ('nsmc_num',), fill_value=np.nan, zlib=False)
        temp.setncatts({'units': units,
                        'long_name': long_name,
                        'missing_value': np.nan,
                        'target': meta.loc[obs, 'Measured'],
                        'nsmc_weight': meta.loc[obs, 'Weight'],
                        'vtype': 'obs'})
        temp[:] = obs_file.loc[obs].values


def _add_convergence(obs_file, nc_file):
    converged = pd.isnull(obs_file.iloc[0]).values
    con = nc_file.createVariable('converged', 'i1', ('nsmc_num',), fill_value=-1, zlib=False)
    con.setncatts({'units': 'boolean',
                   'long_name': 'model_converged',
                   'missing_value': -1,
                   'vtype': 'filter'})
    con[:] = converged


def _add_phis(rec_file, opt_lower_rec, opt_upper_rec, nc_file):
    phis = extractphisummary(rec_file).transpose()
    basenames = ['head', 'vert', 'sfx', 'coast', 'sfo', 'drn', 'total']
    lower_phis = extractoptphi(opt_lower_rec)
    upper_phis = extractoptphi(opt_upper_rec)
    for bn in basenames:
        temp = nc_file.createVariable('phi_{}'.format(bn), 'f8', ('nsmc_num',), fill_value=np.nan, zlib=False)
        temp.setncatts({'units': 'none',
                        'long_name': 'phi for {}'.format(bn),
                        'missing_value': np.nan,
                        'vtype': 'obs'})
        temp_optlower = lower_phis.loc[bn].values
        temp_optupper = upper_phis.loc[bn].values
        temp[:] = np.concatenate((phis.loc[:, bn].values, temp_optlower, temp_optupper))


def _add_filter_1(f1txt, nc_file):
    nsmc_nums = np.array(nc_file.variables['nsmc_num'])
    with open(f1txt) as f:
        pass_nums = [int(e.strip()) for e in f.readlines()] + [-1, -2]
    passed = np.in1d(nsmc_nums, pass_nums).astype(int)
    modfilter = nc_file.createVariable('filter1', 'i1', ('nsmc_num',), fill_value=-1, zlib=False)
    modfilter.setncatts({'units': 'boolean',
                         'long_name': 'filter 1 phi filter',
                         'missing_value': -1,
                         'vtype': 'filter'})
    modfilter[:] = passed


def _add_filter_2(f2txt, nc_file):
    nsmc_nums = np.array(nc_file.variables['nsmc_num'])
    filter1 = np.array(nc_file.variables['filter1'])
    with open(f2txt) as f:
        pass_nums = [int(e.strip()) for e in f.readlines()] + [-1, -2]
    passed = np.in1d(nsmc_nums, pass_nums).astype(int)
    passed[filter1 < 1] = -1

    modfilter = nc_file.createVariable('filter2', 'i2', ('nsmc_num',), fill_value=-1, zlib=False)
    modfilter.setncatts({'units': 'boolean',
                         'long_name': 'filter 2 vert filter',
                         'missing_value': -1,
                         'vtype': 'filter'})
    modfilter[:] = passed


def _add_filter_3(f3txt, nc_file):
    nsmc_nums = np.array(nc_file.variables['nsmc_num'])
    filter1 = np.array(nc_file.variables['filter1'])
    with open(f3txt) as f:
        pass_nums = [int(e.strip()) for e in f.readlines()] + [-1, -2]
    passed = np.in1d(nsmc_nums, pass_nums).astype(int)
    passed[filter1 < 1] = -1

    modfilter = nc_file.createVariable('filter3', 'i2', ('nsmc_num',), fill_value=-1, zlib=False)
    modfilter.setncatts({'units': 'boolean',
                         'long_name': 'filter 3 piezo filter',
                         'missing_value': -1,
                         'vtype': 'filter'})
    modfilter[:] = passed


def _add_filter_4(f4txt, nc_file):
    nsmc_nums = np.array(nc_file.variables['nsmc_num'])
    filter1 = np.array(nc_file.variables['filter1'])
    with open(f4txt) as f:
        pass_nums = [int(e.strip()) for e in f.readlines()] + [-1, -2]
    passed = np.in1d(nsmc_nums, pass_nums).astype(int)
    passed[filter1 < 1] = -1

    modfilter = nc_file.createVariable('filter4', 'i2', ('nsmc_num',), fill_value=-1, zlib=False)
    modfilter.setncatts({'units': 'boolean',
                         'long_name': 'filter 4 intersect of piezo and vert filter',
                         'missing_value': -1,
                         'vtype': 'filter'})
    modfilter[:] = passed


def _add_filter_emma(f5txt, nconv, nc_file):
    emma_data = pd.read_csv(f5txt, index_col=0)
    nsmc_nums = np.array(nc_file.variables['nsmc_num'])
    filters = {}
    filter_atts = {
        'run_mt3d': {'long_name': 'run mt3d',
                     'comments': 'the models which we ran mt3d on'},
        'emma_converge': {'long_name': 'End member mixing converged',
                          'comments': 'there were some random models which did not converge in mt3d denoted as True'},
        'n_converge': {'long_name': 'median nitrate converged',
                       'comments': 'there were some random models which did not converge in mt3d denoted as True'},
        'emma_no_wt': {
            'long_name': 'EMMA with no group weighting individual targets weighted by 1/sd of the EMMA analysis',
            'comments': 'the bottom 10% of data measured by the EMMA phi with no intergroup weights individual targets'
                        ' weighted by 1/sd of the EMMA analysis'},
        'emma_eq_wt': {'long_name': 'EMMA with numerically equally weighting',
                       'comments': 'the bottom 10 of the realisation measured by the emma phi where the inter group '
                                   'weighting was set by 1/ number of obs in each group individual targets weighted by'
                                   ' 1/sd of the EMMA analysis'},
        'emma_chch_wt': {'long_name': 'EMMA with CHCH weighted',
                         'comments': 'bottom 10% with chch groups weighted up by 2 orders of magnitude '
                                     'individual targets weighted by 1/sd of the EMMA analysis'},
        'emma_str_wt': {'long_name': 'EMMA with streams weighted',
                        'comments': 'bottom 10% with stream groups weighted up by 2 orders of magnitude'
                                    ' individual targets weighted by 1/sd of the EMMA analysis'},
        'emma_ewf_wt': {'long_name': 'EMMA with eyrewell forest weighted',
                        'comments': 'bottom 10% with eyrewell forest groups weighted up by 2 orders of magnitude '
                                    'individual targets weighted by 1/sd of the EMMA analysis'},
    }
    # ran mt3d
    ran = np.in1d(nsmc_nums, emma_data.index)
    filters['run_mt3d'] = ran.astype(int)

    # emma converged
    emmaconverged = np.in1d(nsmc_nums, emma_data.loc[emma_data.notnull().any(axis=1)].index)
    f = emmaconverged.astype(int)
    f[~ran] = -1
    filters['emma_converge'] = f

    # N med load converged
    with open(nconv) as f:
        nums = [int(e) for e in f.readlines()]

    f = np.in1d(nsmc_nums, nums).astype(int)
    f[~ran] = -1
    filters['n_converge'] = f

    # emma no weighting
    temp_data = emma_data['no_weighting']
    f = np.in1d(nsmc_nums, emma_data[temp_data <= temp_data.quantile(0.1)].index).astype(int)
    f[~emmaconverged] = -1
    filters['emma_no_wt'] = f

    # emma equal weighting
    temp_data = emma_data['equal_num']
    f = np.in1d(nsmc_nums, emma_data[temp_data <= temp_data.quantile(0.1)].index).astype(int)
    f[~emmaconverged] = -1
    filters['emma_eq_wt'] = f

    # emma chch weighting
    temp_data = emma_data['chch_weighted']
    f = np.in1d(nsmc_nums, emma_data[temp_data <= temp_data.quantile(0.1)].index).astype(int)
    f[~emmaconverged] = -1
    filters['emma_chch_wt'] = f

    # emma stream weighting
    temp_data = emma_data['stream_weighted']
    f = np.in1d(nsmc_nums, emma_data[temp_data <= temp_data.quantile(0.1)].index).astype(int)
    f[~emmaconverged] = -1
    filters['emma_str_wt'] = f

    # emma ewf weighting
    temp_data = emma_data['ewf_weighted']
    f = np.in1d(nsmc_nums, emma_data[temp_data <= temp_data.quantile(0.1)].index).astype(int)
    f[~emmaconverged] = -1
    filters['emma_ewf_wt'] = f

    for key in filters.keys():
        modfilter = nc_file.createVariable(key, 'i2', ('nsmc_num',), fill_value=-1, zlib=False)
        modfilter.setncatts({'units': 'boolean',
                             'long_name': filter_atts[key]['long_name'],
                             'comments': filter_atts[key]['comments'],
                             'missing_value': -1,
                             'vtype': 'filter'})
        modfilter[:] = filters[key]


def make_netcdf_nsmc(nc_outfile, rrffile, rec_file, opt_lower_rec, opt_upper_rec, rch_ppt_tpl, kh_kv_ppt_file,
                     rei_file, opt_lower_rei, opt_upper_rei, pst_file, prior_sds, prior_ksds_dir, postopt_sds,
                     opt_pst_file,
                     f1txt, f2txt, f3txt, f4txt, f5txt, ncontxt, recalc=False):
    """

    :param nc_outfile: path to write the netcdf file
    :param rrffile: the mc rrf file
    :param rec_file: the mc rec file
    :param opt_lower_rec: the lower optimised model rec file (NsmcBase)
    :param opt_upper_rec: the lower optimised model rec file (NsmcBaseB)
    :param rch_ppt_tpl: the recharge pilot point template file (path)
    :param kh_kv_ppt_file: kh or kv pilot point file for layer 1
    :param rei_file: rei file which holds the target values and weights
    :param opt_lower_rei: file which holds the observations for the lower optimised model
    :param opt_upper_rei: file which holds the observation for the upper optimised model
    :param pst_file: the pest contol file for the NSMC
    :param prior_sds: text file with the prior sds for all non khv
    :param prior_ksds_dir: directory with all of the Kh and Kv sds
    :param postopt_sds: textfile with the sds after sensitivity matrix
    :param opt_pst_file: the pest file for the original optimisation to pull kvh priors
    :param f1txt: the text file with a list that passed filter 1 phi
    :param f2txt: the text file with a list that passed filter 2 vert
    :param f3txt: the text file with a list that passed filter 3 piezo
    :param f4txt: the text file with a list that passed filter 4 intersect of 2 & 3
    :param f5txt: the text file with a list that passed filter 5 end member mixing
    :param recalc: if True pull the data out of the rrffile again
    :return:
    """
    # get the data
    obs_file = os.path.join(smt.temp_file_dir, 'temp_obs.csv')
    param_file = os.path.join(smt.temp_file_dir, 'temp_param.csv')
    if not os.path.exists(obs_file) and not recalc:  # to assist with debugging speed
        obs, param = extractrrf(rrffile=rrffile)
        param.to_csv(param_file)
        obs.to_csv(obs_file)
    else:
        obs = pd.read_csv(obs_file, index_col=0)
        param = pd.read_csv(param_file, index_col=0)

    # add optimised model to obs and param, which will then propagate through
    lower_obs = extract_obs_opt_rei(opt_lower_rei, -1)
    upper_obs = extract_obs_opt_rei(opt_upper_rei, -2)

    obs = pd.concat((obs, pd.DataFrame(lower_obs), pd.DataFrame(upper_obs)), axis=1)

    lower_param = param_from_rec(opt_lower_rec, -1)
    upper_param = param_from_rec(opt_upper_rec, -2)
    param = pd.concat((param, pd.DataFrame(lower_param), pd.DataFrame(upper_param)), axis=1)

    pst_param = rd_pst_parameters(pst_file)

    prior_sd_data = pd.read_table(prior_sds, delim_whitespace=True, index_col=0, names=['sd'])
    postopt_sd_data = pd.read_table(postopt_sds, delim_whitespace=True, skiprows=1, index_col=0,
                                    names=['sd'])

    # set up netcdf file
    nc_file = nc.Dataset(nc_outfile, 'w')
    nc_file.notes = """phi lower and phi upper are -1 and -2, respectivly and where present will appear at the end of 
                    the nsmc_variable.  observations with "z" at the end of the name are the weighted 
                    interpretations of zero or near zero targets to make pest notice them (pest interets these).  
                    Their non z brothers are the actual values that were defined as targets
                    all params and obs are in individual variables except: 
                    params: rch pilot points, sfr condances, kh and kv values, and drain conductances
                    obs: head observations (listed as wells)
                    they type of data for each variable is denoted by 'vtype' which can be:
                        'param': parameters for pest
                        'obs': observations
                        'meta': metadata for either parameter or observations (see dimension)
                        'filter': convergence or filters
                        'dim': dimension ids
                        
                    there are two sets of standard deviations for each variable 
                    either as a variable(see params/obs lists above) or attribute :
                        'p_sd': the first standard deviation
                        'j_sd': the standard deviation after the sensitivity matrix (as passed to pest for the NSMC)
                    there is also an attribute (of the base variable or sd variable see above): 
                    'sd_type' either 'lin' or 'log' the true SD or the SD of the log(base 10) transformed data, respectively 
                    the k priors were reduced to 70% for kv and 50% for kh from the true coveriance matrix
                    note that these distributions are not truly bayesian as the distribution is centered by the initial 
                    value and truncated at the bounds
                    
                    'opt_p' is the optimisation prior (e.g. the prior for the optimisation process)
                    """

    # make dimensions
    nc_file.createDimension('nsmc_num', nsmc_dim)
    nc_file.createDimension('layer', layer_dim)
    nc_file.createDimension('rch_ppt', rch_dim)
    nc_file.createDimension('sfr_cond', sfr_dim)
    nc_file.createDimension('khv_ppt', khv_dim)
    nc_file.createDimension('well_name', well_dim)
    nc_file.createDimension('drns', drn_dim)

    # variables
    nsmc_num = nc_file.createVariable('nsmc_num', 'i4', ('nsmc_num',), fill_value=-9, zlib=False)
    nsmc_num.setncatts({'units': 'none',
                        'long_name': 'Null Space Monte Carlo Realisation Number',
                        'comments': 'unique identifier phi lower and phi upper are -1 and -2, respectively',
                        'missing_value': -9,
                        'vtype': 'dim'})
    nsmc_num[:] = range(1, nsmc_dim + -1) + [-1, -2]

    layer = nc_file.createVariable('layer', 'i4', ('layer',), fill_value=-9, zlib=False)
    layer.setncatts({'units': 'none',
                     'long_name': 'model layer',
                     'comments': '1 indexed',
                     'missing_value': -9,
                     'vtype': 'dim'})
    layer[:] = range(1, layer_dim + 1)

    # parameters
    _add_simple_params(param, pst_param, prior_sd_data, postopt_sd_data, nc_file)

    _add_rch_params(param, rch_ppt_tpl, pst_param, prior_sd_data, postopt_sd_data, nc_file)

    _add_sfr_cond(param, pst_param, prior_sd_data, postopt_sd_data, nc_file)

    _add_drain_cond(param, pst_param, prior_sd_data, postopt_sd_data, nc_file)

    _add_kv_kh(param, kh_kv_ppt_file, pst_param, prior_ksds_dir, postopt_sd_data, opt_pst_file, nc_file)

    # add observations convergence and phis
    _add_well_obs(obs, rei_file, nc_file)

    _add_other_obs(obs, rei_file, nc_file)

    _add_convergence(obs, nc_file)

    _add_phis(rec_file, opt_lower_rec, opt_upper_rec, nc_file)

    # add filters
    _add_filter_1(f1txt, nc_file)  # phi filter
    _add_filter_2(f2txt, nc_file)  # vert filter
    _add_filter_3(f3txt, nc_file)  # piezo filter
    _add_filter_4(f4txt, nc_file)  # intersect of vert and piezo filter
    _add_filter_emma(f5txt, ncontxt, nc_file)  # end member mixing filter

    # add general comments
    nc_file.description = (
        'parameters, observations, and filters and metadata for the waimakariri Null Space Monte Carlo.'
        'The nsmc_num are unique ids.  -1 correspons to the lower phi optimisation(NsmcBase) and -2 corresponds to the '
        'upper phi optimisation (NsmcBaseB)')
    nc_file.history = 'created {}'.format(datetime.datetime.now().isoformat())
    nc_file.source = 'script: {}'.format(sys.argv[0])
    nc_file.close()


if __name__ == '__main__':
    data_dir = "{}/from_gns/nsmc".format(smt.sdp)
    make_netcdf_nsmc(nc_outfile=env.gw_met_data("mh_modeling/netcdfs_of_key_modeling_data/nsmc_params_obs_metadata.nc"),
                     rrffile="{}/aw_ex_mc/aw_ex_mc.rrf".format(data_dir),
                     rec_file="{}/aw_ex_mc/aw_ex_mc.rec".format(data_dir),
                     opt_lower_rec="{}/AW_PHILOW_PIEZO/AW_PHILOW_PIEZO/aw_ex_philow_piez.rec".format(data_dir),
                     opt_upper_rec="{}/AW_PHIUPPER_PIEZO/AW_PHIUPPER_PIEZO/aw_ex_phiupper_piez.rec".format(data_dir),
                     rch_ppt_tpl="{}/AW_PHILOW_PIEZO/AW_PHILOW_PIEZO/rch_ppts.tpl".format(data_dir),
                     kh_kv_ppt_file="{}/AW_PHILOW_PIEZO/AW_PHILOW_PIEZO/KV_ppk_01.tpl".format(data_dir),
                     rei_file="{}/aw_ex_jco/aw_ex_jco.rei".format(data_dir),
                     opt_lower_rei="{}/AW_PHILOW_PIEZO/AW_PHILOW_PIEZO/aw_ex_philow_piez.rei".format(data_dir),
                     opt_upper_rei="{}/AW_PHIUPPER_PIEZO/AW_PHIUPPER_PIEZO/aw_ex_phiupper_piez.rei".format(data_dir),
                     pst_file="{}/aw_ex_mc/aw_ex_mc.pst".format(data_dir),
                     prior_sds="{}/sds/PriorSDs.txt".format(data_dir),
                     prior_ksds_dir="{}/sds/ppk_priorSD".format(data_dir),
                     postopt_sds="{}/sds/aw_ex_postopt_sd.txt".format(data_dir),
                     opt_pst_file='{}/from_gns/NsmcBase/AW20171024_2_i2_optver/i2/aw_ex_reg_wtadj_manwtadj_midcal.pst'.format(
                         smt.sdp),
                     f1txt="{}/F1_filtered_parsets.txt".format(data_dir),
                     f2txt="{}/F2_vert_filtered_parsets.txt".format(data_dir),
                     f3txt="{}/F2_piez_filtered_parsets.txt".format(data_dir),
                     f4txt="{}/F2_intersect_filtered_parsets.txt".format(data_dir),
                     f5txt="{}/emma_phis.csv".format(data_dir),
                     ncontxt="{}/mednload_converged.txt".format(data_dir))
    print('done')
