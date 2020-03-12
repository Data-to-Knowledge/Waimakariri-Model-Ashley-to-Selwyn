"""
Author: matth
Date Created: 30/09/2017 9:14 AM
"""

from __future__ import division
from waimak_extended_boundry import smt
import numpy as np
import pickle
import os


def get_zone_array_index(zones, recalc=False): #todo manage data
    """
    returns a boolean array with true for the zones listed
    :param zones: one or more of chch, waimak, selwyn, inland_waimak, coastal_waimak
    :return:
    """
    pickle_path = os.path.join(smt.pickle_dir, 'cwms_zone_arrays.p')

    zones = np.atleast_1d(zones)
    if os.path.exists(pickle_path) and not recalc:
        zones_idx_dict = pickle.load(open(pickle_path))
    else:
        zone_dict = {'waimak': 4, 'chch': 7, 'selwyn': 8}
        new_no_flow = smt.get_no_flow()
        zones_idx = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/cwms_zones.shp".format(smt.sdp), 'ZONE_CODE')
        zones_idx[~new_no_flow[0].astype(bool)] = np.nan
        zones_idx_dict = {}
        for zone in zone_dict.keys():
            use_code = zone_dict[zone]
            temp_idx = zones_idx == use_code
            zones_idx_dict[zone] = temp_idx
        # add waimak coastal/inland
        zones_idx = smt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/coastal_inland_waimak.shp".format(smt.sdp),
                                                  'Id', True)
        temp_idx = (zones_idx==0) & zones_idx_dict['waimak']
        zones_idx_dict['coastal_waimak'] = temp_idx
        temp_idx = (zones_idx==1) & zones_idx_dict['waimak']
        zones_idx_dict['inland_waimak'] = temp_idx
        pickle.dump(zones_idx_dict, open(pickle_path, 'w'))
    out_idx = np.zeros((smt.rows,smt.cols))
    for key in zones:
        out_idx += zones_idx_dict[key]
    out_idx = out_idx.astype(bool)

    return out_idx


if __name__ == '__main__':
    test = get_zone_array_index('waimak')
    print('done')
