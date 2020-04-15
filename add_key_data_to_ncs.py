"""
 Author: Matt Hanson
 Created: 15/04/2020 7:55 AM

 add some odds and ends that are missing datasets to several ncs
 """
import netCDF4 as nc
import flopy_mh as flopy
import numpy as np
from waimak_extended_boundry.extended_boundry_model_tools import smt
from waimak_extended_boundry.model_run_tools.data_extraction.ucn_netcdf import _get_sfr_con_map

if __name__ == '__main__':
    pass

# add comment about dry cells == -888 in heads netcdf
run_hds = False

if run_hds:
    hds_paths = [
        'D:/Waimakariri_model_input_data\\recommended\\post_filter1_hds.nc',
        'D:/Wdc\\data\\Num_GW_model\\netcdfs_of_key_modeling_data\\post_filter1_hds.nc',
        'E:/backup_model_input_data\\Waimakariri_model_input_data_13_4_20\\recommended\\post_filter1_hds.nc',
        'E:/Wdc\\data\\Num_GW_model\\netcdfs_of_key_modeling_data\\post_filter1_hds.nc']

    for p in hds_paths:
        dataset = nc.Dataset(p, 'a')
        dataset.description = '''The heads values for all models which passed filter one (phi filter).
          These were recalculated for filters 2-4
          dry cells are set to -888, see hdry in modflow manual for more details'''

run_add_str_obs = True

if run_add_str_obs:

    mednload_paths = [
        'E:/Wdc\\data\\Num_GW_model\\netcdfs_of_key_modeling_data\\mednload_unc.nc',
        'D:/Wdc\\data\\Num_GW_model\\netcdfs_of_key_modeling_data\\mednload_unc.nc',
        'D:/Waimakariri_model_input_data\\recommended\\post_filter1_mednload_unc.nc',
        'E:/backup_model_input_data\\Waimakariri_model_input_data_13_4_20\\recommended\\post_filter1_mednload_unc.nc'

    ]

    temp_dataset = nc.Dataset(mednload_paths[0])
    nsmc_nums = np.array(temp_dataset.variables['nsmc_num'][:])
    temp_dataset.close()

    nidx = np.where(nsmc_nums == -1)[0][0]
    philow_path = r"C:\matt_modelling_unbackedup\add_to_ncs\mt_aw_ex_mednload_philow_tvd\mt_aw_ex_mednload.ucn"
    philow = flopy.utils.UcnFile(philow_path).get_alldata(nodata=-1)[0]

    # create str obs
    str_obs = np.zeros((len(nsmc_nums), smt.rows, smt.cols)) * np.nan

    for i, n in enumerate(nsmc_nums):
        if i%100 == 0:
            print('reading sobs {} of {}'.format(i, len(nsmc_nums)))
        if n == -1:
            sobpath = r"C:\matt_modelling_unbackedup\add_to_ncs\mt_aw_ex_mednload_philow_tvd\mt_aw_ex_mednload.sobs"
        elif n == -2:
            sobpath = r"C:\matt_modelling_unbackedup\add_to_ncs\MedNload_sobsrepo_files\mt_aw_ex_mednload_phiupper.sobs"
        else:
            sobpath = r"C:\matt_modelling_unbackedup\add_to_ncs\MedNload_sobsrepo_files\mt_aw_ex_mednload_{}.sobs".format(
                n)

        temp = _get_sfr_con_map(sobpath)
        str_obs[i] = temp

    for p in mednload_paths:
        print('writing to {}'.format(p))
        dataset = nc.Dataset(p, 'a')

        # add philow ucn
        dataset.variables['mednload'][nidx] = philow

        # add str_obs
        try:
            dataset.variables['sobs_mednload'][:] = str_obs

        except:
            temp = dataset.createVariable('sobs_mednload', 'f4', ('nsmc_num', 'row', 'col'), fill_value=np.nan,
                                          zlib=False)
            temp.setncatts({'units': 'mg/L',
                            'long_name': 'SFR concentarations for mednload',
                            'missing_value': np.nan})

            temp[:] = str_obs
