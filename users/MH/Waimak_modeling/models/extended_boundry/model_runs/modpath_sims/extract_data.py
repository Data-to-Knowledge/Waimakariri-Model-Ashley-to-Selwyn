# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 6/12/2017 1:37 PM
"""

from __future__ import division
from core import env
import pandas as pd
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.modpath_sims.setup_forward_modpath import \
    part_group_cell_mapper
import netCDF4 as nc
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import itertools
import pickle
import numpy as np
import flopy

# particle id moves to 0 indexed
# particle id links the endpoint file and path file
# make a unique id from a string of k,i,j
# store as a netcdf with 2 variables cellid and percentage
# dimensions of (layer, row, col, ids) +- realisation
# a vectorised dictionary call to extract 'variable' numbers and their corresponding values...
# extraction could take some time, but this is fine as it will be a rare occurance
# the particle
# it looks like every cell occurance that the particle passes through is in teh pathline file maybe confirm with brioch/mike

def open_path_file_as_df(path):
    names = ['Particle_ID',
             'Particle_Group',
             'Time_Point_Index',
             'Cumulative_Time_Step',
             'Tracking_Time',
             'Global_X',
             'Global_Y',
             'Global_Z',
             'Layer',
             'Row',
             'Column',
             'Grid',
             'Local_X',
             'Local_Y',
             'Local_Z',
             'Line_Segment_Index',
             ]
    data = pd.read_table(path, skiprows=3, names=names, delim_whitespace=True)

    # convert to zero indexing
    data.loc[:, 'Layer'] += -1
    data.loc[:, 'Row'] += -1
    data.loc[:, 'Column'] += -1
    return data


def _get_group_num(x):
    return x.iloc[0]


def extract_forward_data(path):
    """
    extract the data and export as pd.DataFrame with index of cell_ref_id, Particle_Group and column of Particle_ID count
    :param path: path to the pathline file
    :return:
    """
    # for now assume that I can hold the full thing in memory, but watch
    drop_names = [
        'Time_Point_Index',
        'Cumulative_Time_Step',
        'Tracking_Time',
        'Global_X',
        'Global_Y',
        'Global_Z',
        'Grid',
        'Local_X',
        'Local_Y',
        'Local_Z',
        'Line_Segment_Index',
    ]
    print('reading data')
    data = open_path_file_as_df(path)
    print('simplifying data')
    data.drop(drop_names, 1, inplace=True)
    # make a ref cell id and make sure it is zero indexed
    data['ref_cell_id'] = ['{:02d}_{:03d}_{:03d}'.format(k, i, j) for k, i, j in
                           data.loc[:, ['Layer', 'Row', 'Column']].itertuples(False, None)]
    data.drop(['Layer', 'Row', 'Column'], axis=1, inplace=True)
    # now for some fancy groupby operations
    print('calculating percentages')
    outdata = data.groupby(['ref_cell_id', 'Particle_ID']).aggregate({'Particle_Group': _get_group_num}).reset_index()
    outdata = outdata.groupby(['ref_cell_id', 'Particle_Group']).count().astype(float)

    # make this output a fraction
    outdata = outdata.rename(columns={'Particle_ID': 'fraction'})
    outdata = outdata.reset_index().set_index('ref_cell_id')

    return outdata


def save_forward_data(path, outpath):
    # save the data extracted above to an emulator netcdf
    # keep the group id to locate cells, but make a linker (e.g. pass the dictionary to the dataframe)
    data = extract_forward_data(path)
    data.to_hdf(outpath, 'emulator', mode='w')


def extract_back_data(path_path, group_mapper_path, hds_path, return_packed_bits=False):
    """

    :param path_path: the pathline file
    :param group_mapper_path: the file to the group mapper produced in set_up_reverse_modpath
    :param hds_path: path to the simulation heads file
    :param return_packed_bits: bool if True return the array as a boolean array and packedbits
    :return:
    """
    # for now assume that I can hold the full thing in memory, but watch
    print('extracting data')
    drop_names = [
        'Time_Point_Index',
        'Cumulative_Time_Step',
        'Tracking_Time',
        'Global_X',
        'Global_Y',
        'Global_Z',
        'Grid',
        'Local_X',
        'Local_Y',
        'Line_Segment_Index',
    ]
    print('reading data')

    # set the local grid maximum so that particles are only counted if they are in the top (1-local_z_max)
    # percent of the cell
    local_z_max = 0.8 #todo this may need to iterativly change scott is critical

    data = open_path_file_as_df(path_path)
    print('simplifying data')
    data.drop(drop_names, 1, inplace=True)
    no_flow = smt.get_no_flow(0)
    hds = flopy.utils.HeadFile(hds_path).get_alldata()[0]
    active_top_cells = smt.model_where((no_flow[np.newaxis] == 1) & (hds[:1] > -777))
    temp = smt.model_where(hds < -777)
    temp = [(e[0]+1, e[1], e[2]) for e in temp]
    active_top_cells.extend(temp)
    # make a ref cell id and make sure it is zero indexed
    data = data.set_index(['Layer', 'Row', 'Column']).loc[active_top_cells]
    data.index.names = ['Layer','Row','Column']
    data.dropna(how='all',inplace=True)
    data.reset_index(inplace=True)  # just keep all data in top active layer
    outdata = {}
    group_mapper = pd.read_csv(group_mapper_path, index_col=0, names=['key', 'val'])['val'].to_dict()
    print("creating maps")
    for i,g in enumerate(set(data.Particle_Group)):
        print('{} of {}'.format(i,len(group_mapper)))
        temp = data.loc[(data.Particle_Group == g) & (data.Local_Z >= local_z_max), ['Row', 'Column']]
        temp = temp.reset_index().groupby(['Row', 'Column']).count().reset_index().values
        temp_out = smt.get_empty_model_grid().astype(int)
        temp_out[temp[:, 0], temp[:, 1]] = temp[:, 2]
        if return_packed_bits:
            temp_out = np.packbits(temp_out>0)
        outdata[group_mapper[g]] = temp_out

    return outdata


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    test = extract_back_data(r"C:\Users\MattH\Desktop\test_reverse_modpath_strong\test_reverse.mppth",
                             r"C:\Users\MattH\Desktop\test_reverse_modpath_strong\test_reverse_group_mapper.csv",
                             "K:\mh_modeling\modflow_dir_for_source\NsmcReal-00001_for_modpath\NsmcReal-00001_for_modpath.hds")
    smt.plt_matrix(test[1] > 0, base_map=True)
    plt.show()
    print('done')
