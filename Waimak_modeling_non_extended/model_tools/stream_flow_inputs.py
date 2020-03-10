"""
Author: matth
Date Created: 22/05/2017 1:35 PM
"""

from __future__ import division
from copy import deepcopy

in_flow_dict = {
    'ashley': 5,
    # ashley northern tribs
    'glentui': 8,
    'garry': 15,
    'bullock': 22,
    'okuku': 24,
    'makerikeri': 28,

    'cust': 7,
    'cust_biwash': 19,

    'eyre': 1,

    'waimak': 3
}

org_flow_vol = { # not used here
    5: 976320.0,  # ashley
    8: 59213.0,  # glentui
    15: 172535.0,  # garry
    22: 79632.0,  # bullock
    24: 547213.0,  # okuku
    28: 74527.0,  # makerikeri
    7: 14688.0,  # cust
    19: 34560.0,  # cust_biwash
    1: 171936.0,  # eyre
    3: 11300000.0,  # waimak
}


def change_inflow(str_spd, changes, set_or_mult):
    """
    change the input data for a stream package
    :param str_spd: original stress period data for one stress period
    :param changes: dictonary of inflow names and values
    :param set_or_mult: set: set to the value in changes; mult: multiply the value in str_spd by teh value in changes
    :return: stress period data (for one stress period)
    """
    outdata = deepcopy(str_spd)
    if set_or_mult == 'set':
        for key in changes.keys():
            seg = in_flow_dict[key]
            outdata['flow'][(outdata['segment'] == seg) & (outdata['reach'] == 1)] = changes[key]
    elif set_or_mult == 'mult':
        for key in changes.keys():
            seg = in_flow_dict[key]
            outdata['flow'][(outdata['segment'] == seg) & (outdata['reach'] == 1)] *= changes[key]
    else:
        raise ValueError('unexpected value for set_or_mult')

    return outdata
