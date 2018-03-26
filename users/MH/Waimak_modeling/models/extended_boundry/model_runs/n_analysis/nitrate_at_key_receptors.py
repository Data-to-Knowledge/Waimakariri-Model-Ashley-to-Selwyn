# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 29/01/2018 9:24 AM
"""

from __future__ import division
from core import env
import pandas as pd
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.con_from_netcdf import \
    calculate_con_from_netcdf_well, calculate_con_from_netcdf_str
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_at_wells import \
    get_con_at_wells
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_from_streams import \
    get_con_at_str, get_samp_points_df
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import \
    get_stocastic_set
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.all_well_layer_col_row import \
    get_all_well_row_col


def get_well_ids():
    """
    the well number and other key info
    :return: pd.Dataframe
    """
    wdc_wells = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\WDC_wells.csv",
        index_col=0)

    all_wells = get_all_well_row_col()

    private_wells = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\PrivateWellZones.csv",
        index_col=0)
    private_wells = pd.merge(private_wells, all_wells.loc[:, ['layer', 'row', 'col']], right_index=True,
                             left_index=True)
    private_wells = private_wells.dropna()
    private_wells = pd.merge(private_wells,
                             pd.DataFrame(all_wells.loc[:, ['depth', 'nztmx', 'nztmy', 'mid_screen_elv']]),
                             how='left', left_index=True, right_index=True)
    idx = private_wells.depth > 50
    private_wells.loc[idx, 'zone_2'] = private_wells.loc[idx, 'Zone_1'] + '_deep'
    private_wells.loc[~idx, 'zone_2'] = private_wells.loc[~idx, 'Zone_1'] + '_shallow'
    out_wells = pd.concat((wdc_wells,private_wells.loc[:,['Zone_1', 'zone_2']]))
    out_wells.loc[:,'private_public'] = 'private'
    out_wells.loc[out_wells.Zone.notnull(),'private_public'] = 'public'

    return out_wells


def get_str_ids():
    """
    get teh stream reaches for the N analysis
    :return: list
    """
    str_sites = ['ashley_sh1',
                 'cust_skewbridge',
                 'cam_bramleys_s',
                 'cam_marshes_s',
                 'courtenay_kaiapoi_s',
                 'kaiapoi_harpers_s',
                 'kaiapoi_island_s',
                 'northbrook_marshes_s',
                 'ohoka_island_s',
                 'saltwater_factory_s',
                 'southbrook_marshes_s',
                 'taranaki_gressons_s',
                 'taranaki_preeces_s',
                 'waikuku_sh1_s',
                 # round 2
                 'ash_ash_est',
                 'ash_est_all',
                 'cam_end_s',
                 'waikuku_end_s',
                 'taranaki_end_s',
                 'saltwater_end_s',
                 'kaiapoi_end']

    return str_sites


# todo make wells the groups of wells rather than the individual? possibly export both sets of data
def get_n_at_points_single_model(outdir, model_id, ucn_file_path, sobs_path, cbc_path, sfo_path):
    """
    saves streams and well data for a transport run from a single model
    :param outdir: directory to save the output in
    :param model_id: the model id
    :param ucn_file_path: the path to the ucn file
    :param sobs_path: the path to the sobs path
    :param cbc_path:  the path to the cbc path
    :param sfo_path: the path to the stream flow file
    :return:
    """
    # run on gw02
    str_sites = get_str_ids()
    str_data = get_con_at_str(sites=str_sites, ucn_file_path=ucn_file_path, sobs_path=sobs_path,
                              cbc_path=cbc_path, sfo_path=sfo_path, kstpkpers=None, rel_kstpkpers=-1)
    str_data.to_csv(os.path.join(outdir, '{}_stream_data.csv'.format(model_id)))
    wells = get_well_ids()
    well_data = get_con_at_wells(well_list=list(set(wells.index)), unc_file_path=ucn_file_path,
                                 kstpkpers=None, rel_kstpkpers=-1, add_loc=True)
    well_data = pd.merge(wells, well_data, left_index=True, right_index=True)
    well_data.to_csv(os.path.join(outdir, '{}_well_data.csv'.format(model_id)))
    zone_sets = [set(well_data.Zone[well_data.Zone.notnull()]),
            set(well_data.Zone_1[well_data.Zone_1.notnull()]),
            set(well_data.zone_2[well_data.zone_2.notnull()])]
    outdata = {}
    for zone_set, key in zip(zone_sets, ['Zone','Zone_1','zone_2']):
        for zone in zone_set:
            idxs = well_data.loc[well_data[key]==zone].index
            temp = well_data.loc[idxs]
            outdata[zone] = temp['con_gm3_kstp0_kper0'].mean()
    outdata = pd.Series(outdata)
    outdata.to_csv(os.path.join(outdir, '{}_grouped_well_data.csv'.format(model_id)))



def get_n_at_points_nc(outdir, nsmc_nums, ucn_var_name='mednload',
                       ucn_nc_path=r"C:\mh_waimak_model_data\mednload_ucn.nc",
                       cbc_nc_path="C:\mh_waimak_model_data\post_filter1_budget.nc",
                       missing_str_obs='raise'):
    """
    pulls out the concentration data from a netcdf file for the given nsmc_nums save both grouped and raw data
    :param outdir: directory to save the data in
    :param nsmc_nums: the nsmc numbers to pull data out for
    :param ucn_var_name: the variable name for nitrate
    :param ucn_nc_path: the path to the ucn netcdf file, must have sobs
    :param cbc_nc_path: the path to the cbc netcdf file
    :param missing_str_obs: action upon missing str_obs, raise Keyerror, warn, or pass
    :return:
    """
    # run on gw02
    str_sites = get_str_ids()
    wells = get_well_ids()

    str_data = calculate_con_from_netcdf_str(nsmc_nums, ucn_nc_path, ucn_var_name, cbc_nc_path, str_sites,
                                             outpath=os.path.join(outdir,'raw_stocastic_set_str_data.csv'),
                                             missing_str_obs=missing_str_obs).describe(
        percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).transpose()
    str_data.to_csv(os.path.join(outdir, 'stocastic_set_strs.csv'))

    all_well_data = calculate_con_from_netcdf_well(nsmc_nums, ucn_nc_path,
                                               ucn_var_name, list(set(wells.index)),
                                               outpath=os.path.join(outdir,'raw_stocastic_set_well_data.csv'))
    well_data = pd.merge(all_well_data.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).transpose(),
                         wells, right_index=True, left_index=True)
    well_data.to_csv(os.path.join(outdir, 'stocastic_set_wells.csv'))
    zone_sets = [set(well_data.Zone[well_data.Zone.notnull()]),
            set(well_data.Zone_1[well_data.Zone_1.notnull()]),
            set(well_data.zone_2[well_data.zone_2.notnull()])]
    outdata = {}
    for zone_set, key in zip(zone_sets, ['Zone','Zone_1','zone_2']):
        for zone in zone_set:
            idxs = well_data.loc[well_data[key]==zone].index
            temp = all_well_data.transpose().loc[idxs].mean()
            outdata[zone] = temp.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    outdata = pd.DataFrame(outdata).transpose()
    outdata.to_csv(os.path.join(outdir, 'stocastic_set_grouped_wells.csv'))


def get_n_ash_opt_stocastic_set(outdir):
    """
    a wrapper the gets teh data from AshOpt and and the stocastic set
    :param outdir:
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    nsmc_nums = get_stocastic_set(False)
    get_n_at_points_nc(outdir, nsmc_nums)
    get_n_at_points_single_model(outdir, model_id='AshOpt',
                                 ucn_file_path=env.gw_met_data(
                                     r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\mt_aw_ex_mednload.ucn"),
                                 sobs_path=env.gw_met_data(
                                     r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\mt_aw_ex_mednload.sobs"),
                                 cbc_path=r"{}\from_gns\AshOpt\AW20180103_Ash0_Opt\AW20180103_Ash0_Opt\mf_aw_ex.cbc".format(
                                     smt.sdp),
                                 sfo_path=r"{}\from_gns\AshOpt\AW20180103_Ash0_Opt\AW20180103_Ash0_Opt\mf_aw_ex.sfo".format(
                                     smt.sdp))


if __name__ == '__main__':
    get_n_ash_opt_stocastic_set(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points"))
