"""
Author: matth
Date Created: 20/06/2017 11:59 AM
"""

from __future__ import division
from core import env
import flopy

def create_dis_package(m):

    dis = flopy.modflow.mfdis.ModflowDis(m,
                                         nlay=17,
                                         nrow=364,
                                         ncol=365,
                                         nper=1,
                                         delr=200,
                                         delc=200,
                                         laycbd=0,
                                         top=_get_top(),
                                         botm=_get_bottom(),
                                         perlen=1,
                                         nstp=1,
                                         tsmult=1,
                                         steady=True,
                                         itmuni=4,
                                         lenuni=2,
                                         unitnumber=719,
                                         xul=1512162.53275,
                                         yul=5215083.5772,
                                         rotation=0.0,
                                         proj4_str=2193)




def _get_top(): #todo
    raise NotImplementedError()

def _get_bottom(): #todo
    raise NotImplementedError