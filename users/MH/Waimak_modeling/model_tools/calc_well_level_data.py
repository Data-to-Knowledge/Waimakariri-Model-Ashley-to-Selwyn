"""
Author: matth
Date Created: 9/06/2017 10:52 AM
"""

from __future__ import division
from core import env
import flopy_mh as flopy
import numpy as np
import pandas as pd


def calc_well_hydrographs(hds_file, well_data, kstpkpers=None):
    """
    extract hydrographic data from a heads file.
    :param hds_file: either the flopy.headfile instance or the path to a heads file
    :param well_data: dataframe of well name col, row, layer.  can include others
    :param kstpkpers: which kstper to do
    :return:
    """

    if isinstance(hds_file,str):
        hds_file = flopy.utils.HeadFile(hds_file)

    if not {'row','col','layer'}.issubset(set(well_data.keys())):
        raise ValueError('expected columns row, col, layer in well_data')

    if kstpkpers is None:
        kstpkpers = hds_file.get_kstpkper()

    kstpkpers = np.atleast_2d(kstpkpers)  # sort out the possiblity of only one kstpkper
    kstpkpers = [tuple(e) for e in kstpkpers]
    # set up multi index pandas dataframe (kstp, kper, columns=all the drain and stream points calculated)

    idx = pd.MultiIndex.from_tuples(kstpkpers, names=['kstp', 'kper'])
    outdata = pd.DataFrame(index=idx, columns=well_data.index.values)
    for kstpkper in kstpkpers:
        kstp = kstpkper[0]
        kper = kstpkper[1]
        # create hydrograph data for drains
        # average the two values in the month as a first stab
        heads = hds_file.get_data(kstpkper=kstpkper)
        heads[np.isclose(heads,-999)] = np.nan
        for well in well_data.index.values:
            layer, row, col = well_data.loc[well,['layer', 'row', 'col']]
            outdata.loc[kstp, kper].loc[well] = heads[layer, row, col]

    return outdata

