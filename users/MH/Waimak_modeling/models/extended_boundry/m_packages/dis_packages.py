"""
Author: matth
Date Created: 20/06/2017 11:59 AM
"""

from __future__ import division
import flopy_mh as flopy

from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt

def create_dis_package(m):
    """
    create and add the dis package
    :param m: a flopy model instance
    :return:
    """
    elv_db = smt.calc_elv_db()
    dis = flopy.modflow.mfdis.ModflowDis(m,
                                         nlay=smt.layers,
                                         nrow=smt.rows,
                                         ncol=smt.cols,
                                         nper=1,
                                         delr=smt.grid_space,
                                         delc=smt.grid_space,
                                         laycbd=0,
                                         top=elv_db[0],
                                         botm=elv_db[1:],
                                         perlen=1,
                                         nstp=1,
                                         tsmult=1,
                                         steady=True,
                                         itmuni=4,
                                         lenuni=2,
                                         unitnumber=719,
                                         xul=smt.ulx,
                                         yul=smt.uly,
                                         rotation=0.0,
                                         proj4_str='EPSG:2193')






