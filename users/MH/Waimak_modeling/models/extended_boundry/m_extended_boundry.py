"""
Author: matth
Date Created: 22/06/2017 10:59 AM
"""

from __future__ import division
from core import env
import flopy
import m_packages
import os
import shutil
from users.MH.Waimak_modeling.supporting_data_path import sdp
from extended_boundry_model_tools import smt

def create_m_extended_boundry(name, dir_path, safe_mode=True, mt3d_link=False, version=smt.model_version, mfv='mfnwt'):
    # sort out paths for the model #todo wrong input
    name = 'm_ex_bd_v{}-{}'.format(version, name)
    dir_path = '{}/m_ex_bd_v{}-{}'.format(os.path.dirname(dir_path), version, os.path.basename(dir_path))
    # remove all previous files in the directory
    if os.path.exists(dir_path):
        if safe_mode:
            cont = input(
                'create_base_modflow_files will delete the directory:\n {} \n continue y/n\n'.format(dir_path)).lower()
            if cont == 'n':
                raise ValueError('script aborted so as not to overwrite {}'.format(dir_path))

        shutil.rmtree(dir_path)
        os.makedirs(dir_path)
    else:
        os.makedirs(dir_path)

    # create the model
    if mfv == 'mf2005':
        ex = "{}/models_exes/mf2005.exe".format(sdp)
    elif mfv == 'mfnwt':
        ex = "{}/models_exes/MODFLOW-NWT_1.1.2/MODFLOW-NWT_1.1.2/bin/MODFLOW-NWT_64.exe".format(sdp)
    else:
        raise ValueError('unexpected modflow version {}'.format(mfv))
    m = flopy.modflow.Modflow(name, version=mfv, exe_name=ex, model_ws=dir_path)

    # add packages
    if version == 'a':
        sfr_version, seg_v, reach_v = smt.sfr_version, smt.seg_v, smt.reach_v
        k_version = smt.k_version
        wel_version = smt.wel_version

    else:
        raise NotImplementedError('model version {} has not yet been defined'.format(version))

    #todo check that these are right
    m_packages.create_dis_package(m)
    m_packages.create_bas_package(m)
    m_packages.create_lay_prop_package(m, k_version)
    m_packages.create_rch_package(m)
    m_packages.create_wel_package(m, wel_version)
    m_packages.create_drn_package(m, wel_version=wel_version, reach_version=reach_v)
    m_packages.create_sfr_package(m, sfr_version, seg_v, reach_v)

    # add simple packages

    # solver
    if mfv =='mfnwt':
        flopy.modflow.mfnwt.ModflowNwt(m,
                                       headtol=0.01,
                                       fluxtol=500,
                                       maxiterout=100,
                                       thickfact=1e-05,
                                       linmeth=1,
                                       iprnwt=0,
                                       ibotav=0,
                                       options='COMPLEX',
                                       Continue=False,
                                       dbdtheta=0.4,
                                       dbdkappa=1e-05,
                                       dbdgamma=0.0,
                                       momfact=0.1,
                                       backflag=1,
                                       maxbackiter=50,
                                       backtol=1.1,
                                       backreduce=0.7,
                                       maxitinner=50,
                                       ilumethod=2,
                                       levfill=5,
                                       stoptol=1e-10,
                                       msdr=15,
                                       iacl=2,
                                       norder=1,
                                       level=5,
                                       north=7,
                                       iredsys=0,
                                       rrctols=0.0,
                                       idroptol=1,
                                       epsrn=0.0001,
                                       hclosexmd=0.0001,
                                       mxiterxmd=50,
                                       unitnumber=714)
    elif mfv == 'mf2005':
        flopy.modflow.mfpcg.ModflowPcg(m, mxiter=50, iter1=1000, npcond=1, hclose=1e-05, rclose=0.1, relax=1.0,
                                       nbpol=0, iprpcg=0, mutpcg=0, damp=1.0, dampt=0.0, ihcofadd=0, extension='pcg',
                                       unitnumber=714)
    else:
        raise ValueError('unexpected modflow version: {}'.format(mfv))

    if mt3d_link:
        mt3d_outunit = 54
        mt3d_outname = '{}_mt3d_link.ftl'.format(m.name)
        link = flopy.modflow.ModflowLmt(m, output_file_name=mt3d_outname, output_file_unit=mt3d_outunit,
                                        package_flows=['sfr'],
                                        unitnumber=21)

    # oc
    flopy.modflow.mfoc.ModflowOc(m, ihedfm=0, iddnfm=0, chedfm=None, cddnfm=None, cboufm=None, compact=True,
                                 stress_period_data={(0, 0): ['save head', 'save drawdown',
                                                              'save budget', 'print budget'],
                                                     (-1, -1): []},
                                 extension=['oc', 'hds', 'ddn', 'cbc', 'ibo'], unitnumber=[22,30,31,32,33])
    #todo handle output in name file might be fixed by re-defining oc



    raise NotImplementedError()
    return m

if __name__ == '__main__':
    create_m_extended_boundry('test',r'C:\Users\MattH\Desktop\test_dir',safe_mode=False)