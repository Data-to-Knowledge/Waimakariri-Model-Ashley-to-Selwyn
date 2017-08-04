"""
Author: matth
Date Created: 20/06/2017 11:59 AM
"""

from __future__ import division
from core import env
import flopy
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt, _get_constant_heads
import numpy as np

def create_bas_package(m):


    bas = flopy.modflow.mfbas.ModflowBas(m,
                                         ibound=smt.get_no_flow(),
                                         strt=create_starting_heads(),
                                         ifrefm=True,
                                         ixsec=False,
                                         ichflg=False,
                                         stoper=None,
                                         hnoflo=-999.99,
                                         unitnumber=1)


def create_starting_heads():
    hds = np.repeat(smt.calc_elv_db()[0][np.newaxis, :, :],
              smt.layers, axis=0)  # set to top of layer 1
    con_heads = _get_constant_heads()
    idx = (np.isfinite(con_heads)) & (~(np.isclose(con_heads,-999)))
    hds[idx] = con_heads[idx]
    if not all(np.isfinite(hds).flatten()):
        raise ValueError('nan values in starting heads')
    return hds

if __name__ == '__main__':
    create_starting_heads()