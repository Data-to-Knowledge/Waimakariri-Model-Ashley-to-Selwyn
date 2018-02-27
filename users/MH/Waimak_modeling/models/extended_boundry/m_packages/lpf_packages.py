"""
Author: matth
Date Created: 20/06/2017 11:58 AM
"""

from __future__ import division
from core import env
import flopy
import numpy as np

def create_lay_prop_package(m, mfv,k_version=1):
    """
    create the layer property package
    :param m: a flopy model instance
    :param mfv: the modflow version mf2005 or mfnwt (e.g. lpf vs upw package)
    :param k_version: the ks to use, this was never really used
    :return:
    """
    layer_type = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    layer_avg = 0
    chani = -1
    layer_vka = 0
    hk = 1 # brioch is pulling these from the SC data (via the pilot points) so will not populate
    hani = 1
    vka = 1 # brioch is pulling these from the SC data (via the pilot points) so will not populate
    ss = 0
    sy = 0

    if mfv == 'mfnwt':
        flopy.modflow.mfupw.ModflowUpw(m,
                                       laytyp=layer_type,
                                       layavg=layer_avg,
                                       chani=chani,
                                       layvka=layer_vka,
                                       laywet=0,
                                       ipakcb=740,
                                       hdry=-888.0,
                                       iphdry=1,
                                       hk=hk,
                                       hani=hani,
                                       vka=vka,
                                       ss=ss,
                                       sy=sy,
                                       vkcb=0.0,
                                       noparcheck=False,
                                       unitnumber=704)
    elif mfv == 'mf2005':
        flopy.modflow.mflpf.ModflowLpf(m,
                                       laytyp=layer_type,
                                       layavg=layer_avg,
                                       chani=chani,
                                       layvka=layer_vka,
                                       laywet=0,
                                       ipakcb=740,
                                       hdry=-888.0,
                                       iwdflg=0,
                                       wetfct=0.1,  # not using
                                       iwetit=1,  # not using
                                       ihdwet=0,  # not using
                                       hk=hk,
                                       hani=hani,
                                       vka=vka,
                                       ss=ss,
                                       sy=sy,
                                       vkcb=0.0,
                                       wetdry=-0.01,  # not using
                                       storagecoefficient=False,
                                       constantcv=False,
                                       thickstrt=False,
                                       nocvcorrection=False,
                                       novfc=False,
                                       unitnumber=704)
    else:
        raise ValueError('unexpected modflow version {}'.format(mfv))

