# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 11/11/2017 10:54 AM
"""

from __future__ import division
import os
from waimak_extended_boundry.nsmc_exploration_results.combine_nsmc_results import *
from waimak_extended_boundry.model_run_tools import get_model_name_path
from future.builtins import input


def make_modflow_netcdfs(hds_nc_path, bud_nc_path, zlib):
    # chekc the postition of phi lower and upper
    if filter2_3_4_num[-2] != -1:
        raise ValueError('expected -1 (phi lower) to second to end')
    if filter2_3_4_num[-1] != -2:
        raise ValueError('expected -2 (phi upper) to be the last')

    path_basenames = []  # this does not contain phi lower or upper
    for i in filter2_3_4_num:
        if i > 0:
            path_basenames.append('mf_aw_ex_{}'.format(i))

    # directory paths
    hds_base_path = env.gw_met_data("mh_modeling/data_from_gns/hdsrepo")
    sfo_base_path = env.gw_met_data("mh_modeling/data_from_gns/sforepo")
    cbc_base_path = env.gw_met_data("mh_modeling/data_from_gns/cbcrepo")

    # phi lower and phi upper hds or pull from get model path
    phi_lower = get_model_name_path('NsmcBase').replace('.nam', '')
    phi_upper = get_model_name_path('NsmcBaseB').replace('.nam', '')

    # flesh out paths and phi lower and upper
    hds_paths = [os.path.join(hds_base_path, '{}.hds'.format(e)) for e in path_basenames] + [phi_lower + '.hds',
                                                                                             phi_upper + '.hds']
    sfo_paths = [os.path.join(sfo_base_path, '{}.sfo'.format(e)) for e in path_basenames] + [phi_lower + '.sfo',
                                                                                             phi_upper + '.sfo']
    cbc_paths = [os.path.join(cbc_base_path, '{}.cbc'.format(e)) for e in path_basenames] + [phi_lower + '.cbc',
                                                                                             phi_upper + '.cbc']

    description = """
    The heads, or cell budget file values for all models which passed filter one (phi filter).  These were recalculated
    for filters 2-4
    """
    # may be commented out as I needed ot re-run one
    # make_hds_netcdf(nsmc_nums=filter2_3_4_num, hds_paths=hds_paths, description=description, nc_path=hds_nc_path, zlib=zlib)
    make_cellbud_netcdf(nsmc_nums=filter2_3_4_num, sfo_paths=sfo_paths, cbc_paths=cbc_paths,
                        description=description, nc_path=bud_nc_path, zlib=zlib)


if __name__ == '__main__':
    cont = input('are you sure you want to re-run make modflow netcdfs it will overwrite and takes some time y/n')
    if cont != 'y':
        raise ValueError('user interuppted process to prevent overwrite')

    # differenes between the two are because I only have so much space on the sever, but teh uncompressed version is so much faster...
    # make_modflow_netcdfs(hds_nc_path=env.gw_met_data("mh_modeling/netcdfs_of_key_modeling_data/post_filter1_hds.nc"),
    #                     bud_nc_path=env.gw_met_data("mh_modeling/netcdfs_of_key_modeling_data/post_filter1_cell_budgets.nc"),
    #                     zlib=True)

    make_modflow_netcdfs(hds_nc_path=r"C:\mh_waimak_model_data\post_filter1_hds.nc",
                         bud_nc_path=r"C:\mh_waimak_model_data\post_filter1_budget.nc",
                         zlib=False)  # on gw02
