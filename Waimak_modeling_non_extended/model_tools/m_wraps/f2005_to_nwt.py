"""
Author: matth
Date Created: 26/05/2017 1:57 PM
"""

from __future__ import division
import flopy_mh as flopy
from users.MH.Waimak_modeling.supporting_data_path import sdp
from base_modflow_wrapper import get_base_mf_ss


def convert_to_nwt(m, iphdry=0):
    """
    convert from 2005 to nwt
    :param m: model to convert
    :param iphdry: iphdry is a flag that indicates whether groundwater head will be set to hdry when the groundwater 
                   head is less than 0.0001 above the cell bottom (units defined by lenuni in the discretization 
                   package). If iphdry=0, then head will not be set to hdry. If iphdry>0, then head will be set to hdry.
                   If the head solution from one simulation will be used as starting heads for a subsequent simulation,
                   or if the Observation Process is used (Harbaugh and others, 2000), then hdry should not be printed
                   to the output file for dry cells (that is, the upw package input variable should be set as iphdry=0).
                   (default is 0)
    :return: 
    """
    # convert LPF to UPW
    lpf = m.lpf
    m.remove_package('LPF')

    m.version = 'mfnwt'
    flopy.modflow.mfupw.ModflowUpw(m,
                                   laytyp=lpf.laytyp,
                                   layavg=lpf.layavg,
                                   chani=lpf.chani,
                                   layvka=lpf.layvka,
                                   laywet=0,
                                   ipakcb=lpf.ipakcb,
                                   hdry=lpf.hdry,
                                   iphdry=iphdry,
                                   hk=lpf.hk,
                                   hani=lpf.hani,
                                   vka=lpf.vka,
                                   ss=lpf.ss,
                                   sy=lpf.sy,
                                   vkcb=lpf.vkcb,
                                   noparcheck=False,
                                   unitnumber=lpf.unit_number[0])
    m.exe_name = "{}/models_exes/MODFLOW-NWT_1.1.2/MODFLOW-NWT_1.1.2/bin/MODFLOW-NWT_64.exe".format(sdp)



if __name__ == '__main__':
    m = get_base_mf_ss()
    n = convert_to_nwt(m)
    print 'print'