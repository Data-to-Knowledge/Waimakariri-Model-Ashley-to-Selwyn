"""
Author: matth
Date Created: 10/04/2017 11:21 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from core.ecan_io import rd_sql, sql_db
import pickle
import os
import flopy_mh as flopy
from warnings import warn
from users.MH.Waimak_modeling.supporting_data_path import sdp, get_org_mod_path

def get_model_x_y (rows=190, cols=365):
    pixelWidth = pixelHeight = 200  # depending how fine you want your raster

    x_min, y_min = 1512162.53275, 5215083.5772 - 200 * rows  # was 5177083.5772

    # mid points of the array
    x = np.linspace(x_min + 100, (x_min + 100 + pixelWidth * (cols - 1)), cols)
    y = np.linspace((y_min + 100 + pixelHeight * (rows - 1)), y_min + 100, rows)

    model_xs, model_ys = np.meshgrid(x, y)
    return model_xs, model_ys

def df_to_array(dataframe, value_to_map, _3d=False, kij = (17, 190, 365)):
    """
    maps a dataframe to a model array
    :param dataframe: dataframe with at least keys (k),i,j and value
    :param value_to_map: key to the value to map onto the array
    :param _3d: boolean if true map to a 3d array with kIJ
    :return: array
    """
    layers,rows,cols = kij
    if _3d:
        if not np.in1d(['k','i','j',value_to_map],dataframe.keys()).all():
            raise ValueError('dataframe missing keys')
        outdata = np.zeros((layers,rows,cols))
        for i in dataframe.index:
            layer, row, col = dataframe.loc[i,['k','i','j']]
            outdata[layer,row,col] = dataframe.loc[i,value_to_map]
    else:
        if not np.in1d(['i','j',value_to_map],dataframe.keys()).all():
            raise ValueError('dataframe missing keys')
        outdata = np.zeros((rows,cols))
        for i in dataframe.index:
            row, col = dataframe.loc[i,['i','j']].astype(int)
            outdata[row,col] = dataframe.loc[i,value_to_map]

    return outdata


def model_where(condition):
    """
    quick function to return list of model indexes from a np.where style condition assumes that the array in the
    condition has the shape (i,j,k) or (j,k)
    :param condition:
    :return:
    """
    idx = np.where(condition)
    idx2 = [i for i in zip(*idx)]
    return idx2


def get_well_postions(well_nums, screen_handling = 'middle', raise_exct = True, error_log_path=None,
                      one_val_per_well = False):
    """
    gets model indexs from a well position
    :param well_nums: value or iterable of well numbers
    :param screen_handling: how to calculate elevation from wells with screen details
                            middle: returns all unique layers for the set
                            all: returns all unique layers intersected by screens(e.g. from top and bottom elevation)
    :param raise_exct: raise exceptions or continue, see error_log_path for handeling if False
    :param error_log_path: path to error log where values are written if raise_exct is True.  if None, exceptions passed
                            as warnings
    :param one_val_per_well: if True, lists of layer values will be condensed and the most frequent layer value will be
                             returned ensuring there is only one elevation per well. Where there is a tie the lowest
                             index layer (e.g. the smallest number) of the tied values will be returned
    :return:
    """
    if error_log_path is not None:
        f = open(error_log_path,'w')

    #check inputs
    if screen_handling not in ['middle','all']:
        raise ValueError('screen_handling not appropriately specified')

    if isinstance(well_nums,str):
        return_value = True
    else:
        return_value = False

    # load supporting information
    # below is a pickled nparray of bottoms (with layer 1 top on the top of the stack
    elv_db = calc_elv_db()
    well_details = rd_sql(**sql_db.wells_db.well_details)
    well_details = well_details.set_index('WELL_NO')
    screen_details = rd_sql(**sql_db.wells_db.screen_details)

    well_nums = np.atleast_1d(well_nums)

    if well_nums.ndim != 1:
        raise ValueError('unexpected input dimensions')

    # initalize output arrays
    layer = np.zeros(well_nums.shape, dtype=object) * np.nan
    row = np.zeros(well_nums.shape, dtype=object) * np.nan
    col = np.zeros(well_nums.shape, dtype=object) * np.nan
    for i, well in enumerate(well_nums):
        try:
            # first get row and column
            ref_level = well_details.loc[well, 'REFERENCE_RL']
            lat = well_details.loc[well, 'NZTMY']
            lon = well_details.loc[well, 'NZTMX']
            ground_level = well_details.loc[well, 'GROUND_RL'] + ref_level
            if pd.isnull(ground_level):
                rt, ct = convert_coords_to_matix(lon, lat)
                ground_level = elv_db[0, rt, ct]  # take as top of model

            screen_num = well_details.loc[well, 'Screens']
            depth = ground_level - well_details.loc[well, 'DEPTH']  # this is from ground level

            # set elevation
            elv = None
            if screen_num == 0:
                elv = depth  # also from ground level
                layer_temp, row_temp, col_temp = convert_coords_to_matix(lon, lat, elv, elv_db)
            elif screen_num == 1:
                top = ground_level - screen_details['TOP_SCREEN'][screen_details['WELL_NO'] == well]
                bottom = ground_level - screen_details['BOTTOM_SCREEN'][screen_details['WELL_NO'] == well]
                if len(top) != 1 or len(bottom) != 1:
                    raise ValueError('well: {} incorrect number of screen values'.format(well))
                if screen_handling == 'middle':
                    elv = (top.iloc[0] + bottom.iloc[0]) / 2  # also from ground level
                    layer_temp, row_temp, col_temp = convert_coords_to_matix(lon, lat, elv, elv_db)
                elif screen_handling == 'all':
                    layer_temp = []
                    lt1, row_temp, col_temp = convert_coords_to_matix(lon, lat, top, elv_db)
                    lt2, row_temp, col_temp = convert_coords_to_matix(lon, lat, bottom, elv_db)
                    layer_temp.extend(range(lt1, lt2 + 1))  # assume the screen is continous
            else:
                layer_temp = []
                for j in range(1, screen_num + 1):
                    top = ground_level - screen_details['TOP_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                       (screen_details['SCREEN_NO'] == j)]
                    bottom = ground_level - screen_details['BOTTOM_SCREEN'][(screen_details['WELL_NO'] == well) &
                                                             (screen_details['SCREEN_NO'] == j)]
                    if len(top) != 1 or len(bottom) != 1:
                        raise ValueError('well: {} incorrect number of screen values for screen {}'.format(well, j))
                    if screen_handling == 'middle':
                        elv = (top.iloc[0] + bottom.iloc[0]) / 2  # also from ground level
                        lt, row_temp, col_temp = convert_coords_to_matix(lon, lat, elv, elv_db)
                        layer_temp.append(lt)
                    elif screen_handling == 'all':
                        lt1, row_temp, col_temp = convert_coords_to_matix(lon, lat, top, elv_db)
                        lt2, row_temp, col_temp = convert_coords_to_matix(lon, lat, bottom, elv_db)
                        layer_temp.extend(range(lt1, lt2+1)) # assume the screen is continous


            if pd.isnull(np.atleast_1d(layer_temp)).sum() != 0:
                raise ValueError('well: {} does not have screen, reference level, '
                                 'ground level, or depth data'.format(well))

            if isinstance(layer_temp,list):
                if one_val_per_well:
                    layer_temp = max(set(layer_temp),key=layer_temp.count)
                else:
                    layer_temp = list(set(layer_temp))
            layer[i] = layer_temp
            row[i] = row_temp
            col[i] = col_temp
        except Exception as val:
            if raise_exct:
                raise Exception('{} {}'.format(well,val))
            elif error_log_path is None:
                warn('{} {}'.format(well,val))
            else:
                f.write('{} {}\n'.format(well,val))

    if raise_exct:
        # check for more missing data
        if pd.isnull(layer).sum() != 0:
            raise ValueError('wells {} do not have data in layer'.format(well_nums[pd.isnull(layer)]))
        elif pd.isnull(row).sum() != 0:
            raise ValueError('wells {} do not have data in row'.format(well_nums[pd.isnull(row)]))
        elif pd.isnull(col).sum() != 0:
            raise ValueError('wells {} do not have data in col'.format(well_nums[pd.isnull(col)]))

    if return_value:
        return layer[0], row[0], col[0]
    else:
        return layer, row, col


def convert_matrix_to_coords(row, col, layer=None, elv_db=None):
    """
    convert from matix indexing to real world coordinates
    :param row:
    :param col:
    :param layer:
    :param elv_db: 3d numpy array of bottom elevation the 0th index is the top layer of the surface of the model
    :return: lon, lat if layer is not specified or lon, lat, elv
    """
    model_xs, model_ys = get_model_x_y()
    if elv_db is None and layer is not None:
        # below is a pickled nparray of bottoms (with layer 1 top on the top of the stack
        elv_db = calc_elv_db()

    col = int(col)
    row = int(row)
    lon = model_xs[0, col]
    lat = model_ys[row, 0]

    if layer is None:
        return lon, lat
    else:
        elv = elv_db[int(layer):int(layer) + 2, int(row), int(col)].mean()
        return lon, lat, elv


def calc_elv_db(recalc=False):
    """
    calculates the elevation database (bottoms with the top of layer one on top or loads pickel
    :param recalc: force recalculates the elv_db from the gns data even if it is not present
    :return: elv_db
    """
    pickle_path =  "{}/inputs/pickeled_files/bottoms.p".format(sdp)
    if (not os.path.exists(pickle_path)) or recalc:
        m = flopy.modflow.Modflow.load('{}.nam'.format(get_org_mod_path())) # using the default model
        elv_db = np.concatenate((m.dis.top.array[np.newaxis,:,:],m.dis.botm.array))
        pickle.dump(elv_db,open(pickle_path,'w'))

    else:
        elv_db = pickle.load(open( pickle_path))

    return elv_db

def convert_coords_to_matix(lon, lat, elv=None, elv_db=None, return_AE=False):
    """
    convert from real world coordinates to model indexes by comparing value to center point of array.  Where the value
    is on an edge defaults to top left corner (e.g. row 1, col 1, layer 1)
    :param lon: longitude NZTM
    :param lat: latitude NZTM
    :param elv: elevation Lyttleton 43 MSL #not sure
    :param elv_db: 3d numpy array of bottom elevation the 0th index is the top layer of the surface of the model
    :param return_AE: if true return (layer,z0), (row, y0), (col, x0) where x0 ect are the preportion along the cell
                      for use with modpath particles
    :return: (row, col) if elv is not presented or (layer, row, col)
    """
    model_xs, model_ys = get_model_x_y()
    if lon > model_xs[0,:].max() or lon < model_xs[0,:].min():
        raise ValueError('coords not in model domain')

    if lat > model_ys[:,0].max() or lat < model_ys[:,:].min():
        raise ValueError('coords not in model domain')

    layer, row, col = None, None, None
    if elv_db is None and elv is not None:
        # below is a pickled np.array of bottoms (with layer 1 top on the top of the stack)
        elv_db = calc_elv_db()
    # find the index of the closest middle point
    # convert lon to col
    col = np.abs(model_xs[0, :] - lon).argmin()

    # convert lat to row
    row = np.abs(model_ys[:, 0] - lat).argmin()

    if return_AE:
        x0 = (100 - (model_xs[0, :] - lon)[col])/200
        y0 = (100 - (model_ys[:, 0] - lat)[row])/200

    if elv is None:
        if return_AE:
            return (row, y0), (col, x0)
        else:
            return row, col
    else:
        # convert elv to layer
        top = elv_db[:-1, row, col]
        bottom = elv_db[1:, row, col]
        if elv > top.max():
            warn('above top of database setting layer to zero')
            layer = 0
            z0 = 0
        else:
            layer = np.where((top > elv) & (bottom <= elv))
            if len(layer[0]) != 1:
                raise ValueError(
                    'returns multiple indexes for layer')
            layer = layer[0][0]
            if return_AE:
                z0 = (top[layer] - elv) / (top[layer]-bottom[layer])
        if return_AE:
            return (layer,z0), (row, y0), (col, x0)
        else:
            return layer, row, col

def calc_z_loc(row, col, elv, elv_db=None):
    if elv_db is None:
        # below is a pickled np.array of bottoms (with layer 1 top on the top of the stack)
        elv_db = calc_elv_db()

    elv = np.atleast_1d(elv)
    out_data = []
    for e in elv:
        top = elv_db[:-1, row, col]
        bottom = elv_db[1:, row, col]
        if e > top.max():
            warn('above top of database setting layer to zero')
            layer = 0
            z0 = 0
        else:
            layer = np.where((top > e) & (bottom <= e))
            if len(layer[0]) != 1:
                raise ValueError(
                    'returns multiple indexes for layer')
            layer = layer[0][0]
            z0 = (top[layer] - e) / (top[layer] - bottom[layer])
        out_data.append((layer, z0))

    return out_data


def get_all_adjacent_cells(cell_locs):

    cell_locs = np.atleast_2d(cell_locs)

    out_locs = []
    for loc in cell_locs:
        out_locs.append(loc)

        k = loc[0]
        i = loc[1]
        j = loc[2]

        if k != 17:
            out_locs.append((k+1,i,j))
        if k != 0:
            out_locs.append((k-1,i,j))
        if i != 190:
            out_locs.append((k,i+1,j))
        if i != 0:
            out_locs.append((k,i-1,j))
        if j != 365:
            out_locs.append((k,i,j+1))
        if j != 0:
            out_locs.append((k,i,j-1))

    return out_locs


if __name__ == '__main__':
    print(get_all_adjacent_cells((1,1,3)))
