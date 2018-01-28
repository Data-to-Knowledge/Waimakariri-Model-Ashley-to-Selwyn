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
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.con_from_netcdf import calculate_con_from_netcdf_well, calculate_con_from_netcdf_str
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_at_wells import get_con_at_wells
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.data_extraction.data_from_streams import get_con_at_str, get_samp_points_df
from single_zone_delination import create_single_zone_indexs
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import get_stocastic_set

def get_well_ids():
    wdc_wells = pd.read_csv(
        r"\\gisdata\projects\SCI\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model build and optimisation\Nitrate\WDC_wells.csv",
        index_col=0)
    return wdc_wells


def get_str_ids():
    str_data = get_samp_points_df()
    str_data = str_data.loc[str_data.m_type == 'source']
    return list(str_data.index)


def get_n_at_points_single_model(outdir, model_id, ucn_file_path, sobs_path, cbc_path, sfo_path):
    # run on gw02
    str_sites = get_str_ids()
    str_data = get_con_at_str(sites=str_sites, ucn_file_path=ucn_file_path, sobs_path=sobs_path,
                              cbc_path=cbc_path, sfo_path=sfo_path, kstpkpers=None, rel_kstpkpers=-1)
    str_data.to_csv(os.path.join(outdir,'{}_stream_data.csv'.format(model_id)))
    wells = get_well_ids()
    well_data = get_con_at_wells(well_list=wells.index, unc_file_path=ucn_file_path,
                                 kstpkpers=None, rel_kstpkpers=-1, add_loc=True)
    out_wells = pd.merge(wells, well_data,left_index=True, right_index=True)
    out_wells.to_csv(os.path.join(outdir,'{}_well_data.csv'.format(model_id)))


def get_n_at_points_nc(outdir, nsmc_nums, ucn_var_name='mednload',
                       ucn_nc_path=r"C:\mh_waimak_model_data\mednload_ucn.nc",
                       cbc_nc_path="C:\mh_waimak_model_data\post_filter1_budget.nc"):
    # run on gw02
    str_sites = get_str_ids()
    wells = get_well_ids()

    str_data = calculate_con_from_netcdf_str(nsmc_nums, ucn_nc_path, ucn_var_name, cbc_nc_path, str_sites,
                                  outpath=None).describe(percentiles=[0.05,0.25,0.5,0.75,0.95]).transpose()
    str_data.to_csv(os.path.join(outdir, 'stocastic_set_strs.csv'))

    well_data = calculate_con_from_netcdf_well(nsmc_nums, ucn_nc_path,
                                               ucn_var_name, wells.index,
                                               outpath=None).describe(percentiles=[0.05,0.25,0.5,0.75,0.95]).transpose()
    well_data = pd.merge(well_data,wells,right_index=True,left_index=True)
    well_data.to_csv(os.path.join(outdir, 'stocastic_set_wells.csv'))

def get_n_ash_opt_stocastic_set(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    nsmc_nums = get_stocastic_set(False)
    get_n_at_points_nc(outdir,nsmc_nums)
    get_n_at_points_single_model(outdir, model_id='AshOpt',
                                 ucn_file_path=env.gw_met_data(r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\mt_aw_ex_mednload.ucn"),
                                 sobs_path=env.gw_met_data(r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\mt_aw_ex_mednload.sobs"),
                                 cbc_path=r"{}\from_gns\AshOpt\AW20180103_Ash0_Opt\AW20180103_Ash0_Opt\mf_aw_ex.cbc".format(smt.sdp),
                                 sfo_path=r"{}\from_gns\AshOpt\AW20180103_Ash0_Opt\AW20180103_Ash0_Opt\mf_aw_ex.sfo".format(smt.sdp))


if __name__ == '__main__':
    get_n_ash_opt_stocastic_set(env.sci(r"Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_results_at_points"))