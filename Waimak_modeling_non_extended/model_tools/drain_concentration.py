"""
Author: matth
Date Created: 26/04/2017 2:15 PM
"""

from __future__ import division
import numpy as np
import pandas as pd
import flopy_mh as flopy
import pickle
import os
import glob
from polygon_to_model_array import shape_file_to_model_array
from users.MH.Waimak_modeling.supporting_data_path import sdp
from warnings import warn

def get_drn_concentration(location, m, con, mt3d_kskper = None,
                          mf_kskper = None, recalc=False):
    """

    :param location: list drain locations
    :param m: mf model
    :param con: either a UCN file object or the path to a UNC file
    :param mt3d_kskper: mt3d time step and time period to include if None uses the last
    :param mf_kskper: mf time step and time period to include if None use the last
    :param recalc: bool recalculate the drain points
    :return: dataframe of concentration and volume for each point
    """
    location = list(np.atleast_1d(location))
    outdata = pd.DataFrame(index=location, columns = ['con', 'vol'])
    drn_points = get_drn_samp_pts_dict(recalc=recalc)
    if not set(location).issubset(drn_points.keys()):
        raise ValueError('{} has not been added to drain site dictionary'.format(set(location)-set(drn_points.keys())))

    # load cell by cell and concentration data
    cbc = flopy.utils.CellBudgetFile('{}/{}'.format(m.model_ws,m.get_output(unit=741)
                                                              ))
    if mf_kskper is None:
        mf_kskper = cbc.get_kstpkper()[-1] # use the last timestep/time period
    drn_recharge_data = cbc.get_data(kstpkper=mf_kskper,full3D=True)[0][0] #keep only top layer
    if isinstance(con, str):
        con = flopy.utils.UcnFile(con)
    if mt3d_kskper is None:
        mt3d_kskper = con.get_kstpkper()[-1] # use the last time step is none
    gw_conc_data = con.get_data(mt3d_kskper)[0]
    # set null values to nan
    gw_conc_data[np.isclose(gw_conc_data, 1e+30)] = np.nan

    for loc in location:
        mask_array = drn_points[loc]

        if drn_recharge_data[mask_array].mask.sum() != 0:
            warn(('masked values returned for {}'.format(loc)))

        water_volume = drn_recharge_data[mask_array].data.sum()
        load = (gw_conc_data[mask_array] * drn_recharge_data[mask_array].data).sum()
        outdata.loc[loc,'con'] = load/water_volume
        outdata.loc[loc,'vol'] = water_volume

    return outdata


def get_drn_samp_pts_dict(recalc=False):
    """
    gets a dictionary of boolean arrays for each sampling point.  These are ultimately derived from shape files, but
    if possible this function will load a pickled dictionary
    :param recalc: bool if True then the pickled dictionary (if any) will not be re-loaded and instead the dictionary
                   will be calculated from all avalible shapefiles (in base_shp_path)
    :return: dictionary {location: masking array}
    """
    pickle_path = "{}/inputs/pickeled_files/drn_samp_pts.p".format(sdp)
    if os.path.exists(pickle_path) and not recalc:
        drn_con_samp_pts = pickle.load(open(pickle_path))
        return drn_con_samp_pts

    # load all shapefiles in base_shp_path
    base_shp_path = "{}/inputs/shp_files/drain_catchments2/*.shp".format(sdp)
    temp_lst = glob.glob(base_shp_path)
    temp_kys = [os.path.basename(e).split('.')[0] for e in temp_lst]

    shp_dict = dict(zip(temp_kys,temp_lst))

    drn_con_samp_pts = {}
    for loc in shp_dict.keys():
        temp = shape_file_to_model_array(shp_dict[loc],'elv',alltouched=True)
        temp2 = np.zeros(temp.shape,dtype=bool)
        temp2[np.isfinite(temp)] = True
        temp2[np.isnan(temp)] = False
        drn_con_samp_pts[loc] = temp2

    pickle.dump(drn_con_samp_pts,open(pickle_path,mode='w'))
    return drn_con_samp_pts

def get_drn_samp_pts():
    keys = get_drn_samp_pts_dict().keys()
    for key in keys:
        if 'str_' in key:
            keys.remove(key)

    keys.remove('all_drains')
    keys.remove('waimak_drn')
    return keys

if __name__ == '__main__':
    get_drn_samp_pts_dict(True)