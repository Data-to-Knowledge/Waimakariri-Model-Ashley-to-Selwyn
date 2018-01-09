# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/01/2018 3:36 PM
"""

from __future__ import division
from core import env
import flopy
import os

# todo is from the old model tools need to refine it
#todo set up to pass without a modflow model just an FTL
def create_mt3d(m, mt3d_name=None,
                ssm_crch=None, ssm_stress_period_data=None,
                adv_sov=0, adv_percel=1,
                btn_porsty=0.05, btn_scon=0, btn_nprs=0, btn_timprs=None,
                dsp_lon=0.1, dsp_trpt=0.1, dsp_trpv=0.01,
                nper=None, perlen=None, nstp=None, tsmult=None,  # these can be either value of list of values
                                                                             # and must match flow model if it is not SS
                ssflag=None, dt0=0, mxstrn=50000, ttsmult=1.0, ttsmax=0, # these can be either value of list of values
                gcg_isolve = 1, gcg_inner=50, gcg_outer=1):
    """

    :param m: modflow model or path to namefile of a model which has been run (then the model is loaded)
    :param ssm_crch: the recharge concentration for species 1
    :param ssm_stress_period_data: stress period data for other source/sinks
    :param mt3d_name: the name for all mt3d files if none name will mirror that of the modflow model name
    :param adv_sov: is an integer flag for the advection solution option. MIXELM = 0, the standard finite-difference
                    method with upstream or central-in-space weighting, depending on the value of NADVFD; = 1,
                    the forward-tracking method of characteristics (MOC); = 2, the backward-tracking modified method
                    of characteristics (MMOC); = 3, the hybrid method of characteristics (HMOC) with MOC or MMOC
                    automatically and dynamically selected; = -1, the third-order TVD scheme (ULTIMATE).
                    set up for -1 or 0
    :param adv_percel: PERCEL is the Courant number (i.e., the number of cells, or a fraction of a cell)
                       advection will be allowed in any direction in one transport step. For implicit finite-difference
                       or particle-tracking-based schemes, there is no limit on PERCEL, but for accuracy reasons, it is
                       generally not set much greater than one. Note, however, that the PERCEL limit is checked over the
                       entire model grid. Thus, even if PERCEL > 1, advection may not be more than one cell's length at
                       most model locations. For the explicit finite-difference or the third-order TVD scheme,
                       PERCEL is also a stability constraint which must not exceed one and will be automatically reset
                       to one if a value greater than one is specified.
    :param btn_porsty: porosity for the model
    :param btn_scon: float, array of (nlay, nrow, ncol), or filename, or a list (length ncomp) of these for
                     multi-species simulations The starting concentration for the solute transport simulation
    :param btn_nprs: A flag indicating (i) the frequency of the output and (ii) whether the output frequency is
                     specified in terms of total elapsed simulation time or the transport step number. If nprs > 0
                     results will be saved at the times as specified in timprs; if nprs = 0, results will not be saved
                     except at the end of simulation; if NPRS < 0, simulation results will be saved whenever the number
                     of transport steps is an even multiple of nprs. (default is 0).
    :param btn_timprs: The total elapsed time at which the simulation results are saved. The number of entries in timprs
                       must equal nprs. (default is None).
    :param dsp_lon: the longitudinal dispersivity, for every cell of the model grid (unit, L).
    :param dsp_trpt: is a 1D real array defining the ratio of the horizontal transverse dispersivity to the longitudinal
                     dispersivity. Each value in the array corresponds to one model layer.
    :param dsp_trpv: is the ratio of the vertical transverse dispersivity to the longitudinal dispersivity. Each value
                     in the array corresponds to one model layer. Some recent field studies suggest that TRPT is
                     generally not greater than 0.01. Set TRPV equal to TRPT to use the standard isotropic dispersion
                     model (Equation 10 in Chapter 2). Otherwise, the modified isotropic dispersion model is used
    :param ssflag:  If SSFlag is set to SSTATE (case insensitive), the steady-state transport simulation is
                    automatically activated. (see mt3dms_V5_ supplemental for more info) must be an iterable otherwise
                    only the first letter will be written

    :return: mt3d instance
    """
    # create a general MT3d class instance to run most of our mt3d models

    # load model if m is a string
    if isinstance(m,str):
        m = flopy.modflow.Modflow.load(m, model_ws=os.path.dirname(m),forgive=False,
                                       exe_name="{}/models_exes/mf2005.exe".format(sdp))
        if not os.path.exists('{}/{}'.format(m.model_ws,m.output_fnames[0])):
            raise ValueError('Model must have been run prior to mt3d load')
        else:
            print('model has loaded successfully')
    elif isinstance(m, flopy.modflow.mf.Modflow):
        pass
    else:
        raise ValueError('unexpected input for m')


    # packages I'll likely need
    if mt3d_name is None:
        mt3d_name = '{}_mt3d'.format(m.name)
    mt3d = flopy.mt3d.Mt3dms(modelname=mt3d_name,
                             modflowmodel=m,
                             ftlfilename=m.lmt6.output_file_name, #todo could define here
                             ftlfree=False,  # default but right I think
                             version='mt3d-usgs',
                             exe_name="{}/models_exes/mt3d-usgs_Distribution/bin/MT3D-USGS_64.exe".format(sdp), #todo update with brioch's compilation
                             structured=True,  # defualt
                             listunit=500,
                             ftlunit=501,
                             model_ws=m.model_ws, #todo perhaps change
                             load=True,  # defualt
                             silent=0  # defualt
                             )

    # ADV
    if adv_sov >= 1:
        raise ValueError('mt3d object not configured for specified sover {}'.format(adv_sov))

    adv = flopy.mt3d.Mt3dAdv(mt3d,
                             mixelm=adv_sov,
                             percel=adv_percel,
                             mxpart=800000,  # not using particles
                             nadvfd=1,  # default to upstream weighting
                             itrack=3,  # not using particles
                             wd=0.5,  # not using particles
                             dceps=1e-05,  # defualt
                             nplane=2,  # not using particles
                             npl=10,  # not using particles
                             nph=40,  # not using particles
                             npmin=5,  # not using particles
                             npmax=80,  # not using particles
                             nlsink=0,  # not using particles
                             npsink=15,  # not using particles
                             dchmoc=0.0001,  # not using MOC or MMOC or HMOC
                             unitnumber=502
                             )

    # BTN
    #    'BTN error. Required input is None, but no modflow model.'\
    #    ' If no modflow model is passed to Mt3dms, then values '' \
    #    ''must be specified in the BTN constructor for: '' \
    #    ''nlay, nrow, ncol, nper, laycon, delr, delc, htop, dz, '' \
    #    ''perlen, nstp, and tsmult.'

    btn = flopy.mt3d.Mt3dBtn(mt3d,
                             MFStyleArr=False,  # defualt it's a reader, should hopefully not cause problems
                             DRYCell=True,  # pass through dry cells
                             Legacy99Stor=False,  # defualt
                             FTLPrint=False,  # defualt
                             NoWetDryPrint=False,  # defualt shouldn't be a problem
                             OmitDryBud=True,  # as passing through dry cells
                             AltWTSorb=False,  # defualt not using sorbing to my knowledge
                             ncomp=1,  # number of species
                             mcomp=1,  # number of moblile species
                             tunit='D',
                             lunit='M',
                             munit='G',
                             prsity=btn_porsty,
                             icbund=1,  # all cells active
                             sconc=btn_scon,
                             cinact=-9999999,  # defualt
                             thkmin=0.01,  # defualt

                             # printing flags 0 is not print #todo I could probably get away with not printing anything
                             ifmtcn=1,
                             ifmtnp=1,
                             ifmtrf=1,
                             ifmtdp=1,

                             savucn=True,  # default
                             nprs=btn_nprs,
                             timprs=btn_timprs,
                             obs=None,  # default not using as it is easier to pull from the UCN file
                             nprobs=1,  # not using obs so doesnt matter
                             chkmas=True,
                             nprmas=1,  # defualt print mass balance for each time period
                             nper=nper,
                             perlen=perlen,
                             nstp=nstp,
                             tsmult=tsmult,
                             ssflag=ssflag,
                             dt0=dt0,
                             mxstrn=mxstrn,
                             ttsmult=ttsmult,
                             ttsmax=ttsmax,
                             species_names=['N'],
                             extension='btn',
                             unitnumber=503
                             )

    # DSP
    dsp = flopy.mt3d.Mt3dDsp(mt3d,
                             al=dsp_lon,
                             trpt=dsp_trpt,
                             trpv=dsp_trpv,
                             dmcoef=1e-09,  # default don't think I need as only if multidiff True
                             extension='dsp',
                             multiDiff=False,  # only one component
                             unitnumber=504)

    # SSM
    #warnings.warn('SSM Package: mxss is None and modflowmodel is ' +
    #              'None.  Cannot calculate max number of sources ' +
    #              'and sinks.  Estimating from stress_period_data. ')

    ssm = flopy.mt3d.Mt3dSsm(mt3d,
                             crch=ssm_crch,
                             cevt=None,
                             mxss=None,  # default max number of sources and sinks this is calculated from modflow model
                             stress_period_data=ssm_stress_period_data,
                             dtype=None,  # default I should not need to specify this, but we'll see
                             extension='ssm',
                             unitnumber=505,
                             )
    kwargs = {} #todo place holder
    mtsft = flopy.mt3d.Mt3dSft(mt3d,
                               nsfinit=0,
                               mxsfbc=0,
                               icbcsf=0,
                               ioutobs=None,
                               ietsfr=0,
                               isfsolv=1,
                               wimp=0.50,
                               wups=1.00,
                               cclosesf=1.0E-6,
                               mxitersf=10,
                               crntsf=1.0,
                               iprtxmd=0,
                               coldsf=0.0,
                               dispsf=0.0,
                               nobssf=0,
                               obs_sf=None,
                               sf_stress_period_data=None,
                               unitnumber=None,
                               filenames=None,
                               dtype=None,
                               extension='sft',
                               **kwargs)
    #todo add mtsft


    # GCG
    gcg = flopy.mt3d.Mt3dGcg(mt3d,
                             mxiter=gcg_outer,  # defualt
                             iter1=gcg_inner,  # defualt
                             isolve=gcg_isolve,  # 1, Jacobi = 2, SSOR = 3, Modified Incomplete Cholesky (MIC)
                             ncrs=0,  # lump despersion tensor to RHS
                             accl=1,  # defualt and likely not used
                             cclose=1e-05,  # defualt
                             iprgcg=0,  # defualt print max changes at end of each iteration
                             extension='gcg',
                             unitnumber=506,
                             )

    return mt3d


default_mt3d_dict = {
    'm': None, 'mt3d_name': None,
    'ssm_crch': None, 'ssm_stress_period_data': None,
    'adv_sov': 0, 'adv_percel': 1,
    'btn_porsty': 0.05, 'btn_scon': 0, 'btn_nprs': 0, 'btn_timprs': None,
    'dsp_lon': 0.1, 'dsp_trpt': 0.1, 'dsp_trpv': 0.01,
    'nper': None, 'perlen': None, 'nstp': None, 'tsmult': None, 'ssflag': None,
    'dt0': 0, 'mxstrn': 50000, 'ttsmult': 1.0, 'ttsmax': 0
}
