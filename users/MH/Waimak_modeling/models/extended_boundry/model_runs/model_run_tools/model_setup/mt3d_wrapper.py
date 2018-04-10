# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/01/2018 3:36 PM
"""

from __future__ import division
from core import env
import flopy
import os
import shutil
import numpy as np
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_bc_data.wells import \
    get_race_data
import pandas as pd


def create_mt3d(ftl_path, mt3d_name, mt3d_ws,
                ssm_crch=None, ssm_stress_period_data=None, sft_spd=None,
                adv_sov=0, adv_percel=1,
                btn_porsty=0.05, btn_scon=0, btn_nprs=0, btn_timprs=None,
                dsp_lon=0.1, dsp_trpt=0.1, dsp_trpv=0.01,
                nper=1, perlen=1, nstp=1, tsmult=1,  # these can be either value of list of values
                # and must match flow model if it is not SS
                ssflag=None, dt0=0, mxstrn=50000, ttsmult=1.0, ttsmax=0,  # these can be either value of list of values
                gcg_isolve=1, gcg_inner=50, gcg_outer=1):
    """

    :param ftl_path: path to the FTL file to use with MT3D
    :param mt3d_name: the name for all mt3d files if none name will mirror that of the modflow model name
    :param mt3d_ws: working directory for the MT3D model
    :param ssm_crch: the recharge concentration for species 1
    :param ssm_stress_period_data: stress period data for other source/sinks
    :param sft_spd: sft stress period data

    below here can generally be added by the default mt3d dictionary
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

    # check that FTL is in the model_ws folder and if not move it there
    if not os.path.exists(mt3d_ws):
        os.makedirs(mt3d_ws)

    ftl_name = os.path.basename(ftl_path)
    if not os.path.dirname(ftl_path) == mt3d_ws:
        shutil.copyfile(ftl_path, os.path.join(mt3d_ws, ftl_name))

    # packages I'll likely need
    mt3d = flopy.mt3d.Mt3dms(modelname=mt3d_name,
                             modflowmodel=None,
                             ftlfilename=ftl_name,
                             ftlfree=True,  # formatted FTL to handle bug
                             version='mt3d-usgs',
                             exe_name="{}/models_exes/mt3d_usgs_brioch_comp/"  # standard compilation did not converge
                                      "mt3d-usgs-1.0.exe".format(os.path.dirname(smt.sdp)),
                             structured=True,
                             # defualt probably fine, though a point of weakness I don't know what it is
                             listunit=500,
                             ftlunit=501,
                             model_ws=mt3d_ws,
                             load=True,  # defualt
                             silent=0  # defualt
                             )

    # BTN
    elv_db = smt.calc_elv_db()

    btn = flopy.mt3d.Mt3dBtn(mt3d,  # todo btn should be good
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

                             # printing flags 0 is not print
                             ifmtcn=0,
                             ifmtnp=0,
                             ifmtrf=0,
                             ifmtdp=0,

                             savucn=True,  # default
                             nprs=btn_nprs,
                             timprs=btn_timprs,
                             obs=None,  # default not using as it is easier to pull from the UCN file
                             nprobs=1,  # not using obs so doesnt matter
                             chkmas=True,
                             nprmas=1,  # defualt print mass balance for each time period

                             # modflow model parameters
                             nper=nper,
                             perlen=perlen,
                             nstp=nstp,
                             tsmult=tsmult,
                             ncol=smt.cols,
                             nlay=smt.layers,
                             nrow=smt.rows,
                             laycon=[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             delr=smt.grid_space,
                             delc=smt.grid_space,
                             htop=elv_db[0],
                             dz=elv_db[0:-1] - elv_db[1:],

                             ssflag=ssflag,
                             dt0=dt0,
                             mxstrn=mxstrn,
                             ttsmult=ttsmult,
                             ttsmax=ttsmax,
                             species_names=['N'],
                             extension='btn',
                             unitnumber=503
                             )

    # ADV
    if adv_sov >= 1:
        raise ValueError('mt3d object not configured for specified advection solver {}'.format(adv_sov))

    adv = flopy.mt3d.Mt3dAdv(mt3d,  # todo done
                             mixelm=adv_sov,
                             percel=adv_percel,
                             mxpart=5000,  # not using particles
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
    # DSP
    dsp = flopy.mt3d.Mt3dDsp(mt3d,  # todo done
                             al=dsp_lon,
                             trpt=dsp_trpt,
                             trpv=dsp_trpv,
                             dmcoef=0,
                             # default don't think I need as only if multidiff True
                             extension='dsp',
                             multiDiff=False,  # only one component
                             unitnumber=504)

    # SSM
    # warnings.warn('SSM Package: mxss is None and modflowmodel is ' +
    #              'None.  Cannot calculate max number of sources ' +
    #              'and sinks.  Estimating from stress_period_data. ')

    ssm = flopy.mt3d.Mt3dSsm(mt3d,  # todo done
                             crch=ssm_crch,
                             cevt=None,
                             mxss=100000,
                             # default max number of sources and sinks this is calculated from modflow model # define aprioi
                             stress_period_data=ssm_stress_period_data,
                             dtype=None,  # default I should not need to specify this, but we'll see
                             extension='ssm',
                             unitnumber=505,
                             )

    # MTSFT
    mtsft = flopy.mt3d.Mt3dSft(mt3d,  # todo should be done
                               nsfinit=1228,
                               mxsfbc=1228,
                               icbcsf=0,
                               ioutobs=81,
                               ietsfr=0,
                               isfsolv=1,
                               wimp=0.50,
                               wups=1.00,
                               cclosesf=1.0E-6,
                               mxitersf=10,
                               crntsf=1.0,
                               iprtxmd=0,
                               coldsf=0.0,  # inital concentrations in the stream network... should be fine at 0
                               dispsf=0.0,
                               nobssf=1228,
                               obs_sf=range(1, 1228 + 1),  # I think this is the reach set up
                               sf_stress_period_data=sft_spd,
                               unitnumber=None,
                               filenames=None,
                               dtype=None,
                               extension='sft')

    # GCG
    gcg = flopy.mt3d.Mt3dGcg(mt3d,
                             mxiter=gcg_outer,  # defualt
                             iter1=gcg_inner,  # defualt
                             isolve=gcg_isolve,  # 1, Jacobi = 2, SSOR = 3, Modified Incomplete Cholesky (MIC)
                             ncrs=0,  # lump despersion tensor to RHS
                             accl=1,  # defualt and likely not used
                             cclose=1e-06,  # defualt
                             iprgcg=0,  # defualt print max changes at end of each iteration
                             extension='gcg',
                             unitnumber=506,
                             )

    return mt3d


def get_ssm_stress_period_data(wil_race_con=0.1, upper_n_bnd_flux_con=0.1, lower_n_bnd_flux_con=2,
                               well_stream_con=0, llrzf_con=0, ulrzf_con=0, s_race_con=0, chb_con=0):
    """
    get data for teh ssm package from wells.  cmp pathway is the default values
    :param wil_race_con: concentration to apply to will race
    :param upper_n_bnd_flux_con: from the gorge to ~ glentui terrace
    :param lower_n_bnd_flux_con: from the glen tui terrace to broad road
    :param well_stream_con: concentration to apply for the southern streams that are represented by injection wells
    :param llrzf_con: the lower southwestern boundary flux concentration (below sh1)
    :param ulrzf: the upper southwestern boundary flux concentration (above sh1)
    :param s_race_con: the southern water races concentration
    :param chb_con: the concentration for the constant head boundaries
    :return:
    """
    cbc_cells = np.array(smt.model_where(smt.get_no_flow() < -0.1))
    all_wells = get_race_data('NsmcBase')  # the model id parameter does not matter here as I only need the locations
    all_wells.loc[all_wells.nsmc_type == 'llrzf', 'css'] = llrzf_con
    all_wells.loc[all_wells.nsmc_type == 'ulrzf', 'css'] = ulrzf_con
    all_wells.loc[all_wells.nsmc_type == 's_race', 'css'] = s_race_con
    all_wells.loc[all_wells.nsmc_type == 'sriv', 'css'] = well_stream_con
    all_wells.loc[all_wells.nsmc_type == 'n_race', 'css'] = wil_race_con
    all_wells.loc[(all_wells.nsmc_type == 'nbndf') & (all_wells.col <= 150), 'css'] = upper_n_bnd_flux_con
    all_wells.loc[(all_wells.nsmc_type == 'nbndf') & (all_wells.col > 150), 'css'] = lower_n_bnd_flux_con

    num_recs = len(all_wells) + len(cbc_cells)
    ssm_spd = np.recarray((num_recs), flopy.mt3d.Mt3dSsm.get_default_dtype())
    k = np.concatenate((cbc_cells[:, 0], all_wells.loc[:, 'layer'].values))
    i = np.concatenate((cbc_cells[:, 1], all_wells.loc[:, 'row'].values))
    j = np.concatenate((cbc_cells[:, 2], all_wells.loc[:, 'col'].values))
    itype = np.concatenate((np.repeat([1], len(cbc_cells)), np.repeat([2], all_wells.shape[0])))
    css = np.concatenate((np.repeat([chb_con], len(cbc_cells)), all_wells.loc[:, 'css'].values))
    for nm in ['k', 'i', 'j', 'itype', 'css']:
        ssm_spd[nm] = eval(nm)
    ssm_spd = pd.DataFrame(ssm_spd)

    # remove all zero concentrations (should be fine)
    ssm_spd = ssm_spd.loc[~np.isclose(ssm_spd.css, 0)] #todo confirm you do not need to specify the CHB cells
    return ssm_spd.to_records(index=False)


def get_sft_stress_period_data(eyre=0.35, waimak=0.1, ash_gorge=0.1, cust=0.35, cust_biwash=0.1, ash_tribs=0.35,
                               glen_tui=None, garry=None, bullock=None, okuku=None, makerikeri=None):
    ashtrib_input = ash_tribs is not None and all([e is None for e in [glen_tui, garry, bullock, okuku, makerikeri]])
    individual_trib_input = ash_tribs is None and all(
        [e is not None for e in [glen_tui, garry, bullock, okuku, makerikeri]])
    assert ashtrib_input or individual_trib_input, 'ashley tribs must be either set by ash_tribs or individually, not done right'
    if ash_tribs is not None:
        glen_tui = ash_tribs
        garry = ash_tribs
        bullock = ash_tribs
        okuku = ash_tribs
        makerikeri = ash_tribs

    sft_spd = np.recarray((10), flopy.mt3d.Mt3dSft.get_default_dtype())
    nodes = [
        # main inflows
        1,  # eyre inflow
        77,  # waimakariri inflow
        262,  # ashley inflow
        347,  # cust inflow
        # ashley tribs

        373,  # glen tui
        717,  # garry
        822,  # cust biwash
        867,  # bullock creek
        896,  # okuku river
        925,  # makerikeri river

    ]
    scon_dict = {
        1: eyre,
        77: waimak,
        262: ash_gorge,
        347: cust,
        373: glen_tui,
        717: garry,
        822: cust_biwash,
        867: bullock,
        896: okuku,
        925: makerikeri,
    }
    scons = [scon_dict[e] for e in nodes]
    sft_spd['node'] = nodes  # reach ids
    sft_spd['isfbctyp'] = 0  # set all to headwaters sites
    sft_spd['cbcsf0'] = scons

    return sft_spd


def get_default_mt3d_kwargs():
    default_mt3d_dict = {
        'adv_sov': 0,  # matches brioch
        'adv_percel': 1,  # matcheds brioch
        'btn_porsty': 1,  # modified to match brioch
        'btn_scon': 0.1,  # modified to match brioch
        'btn_nprs': 0,  # output timing (only at end)
        'btn_timprs': None,  # not needed
        'dsp_lon': 10,  # modified to match brioch
        'dsp_trpt': 0.1,  # modified to match brioch
        'dsp_trpv': 0.01,  # modified to match brioch
        'nper': 1,  # from modflow model
        'perlen': 7.3050E5,  # modified to match brioch's
        'nstp': 1,  # modified to match briochs
        'tsmult': 1,  # modified to match briochs
        'ssflag': None,  # DO NOT SET
        'dt0': 1,  # modified to match briochs
        'mxstrn': 1000000,  # modified to match briochs
        'ttsmult': 1.1,  # modified to match briochs
        'ttsmax': 50,  # modified to match briochs
        'gcg_isolve': 3,  # modified to match briochs
        'gcg_inner': 500,  # modified to match briochs
        'gcg_outer': 100  # modified to match briochs
    }
    return default_mt3d_dict


if __name__ == '__main__':
    from users.MH.Waimak_modeling.models.extended_boundry.model_runs.model_run_tools.model_setup.realisation_id import get_model
    rch_path = r"K:\mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\nconc_cmp_200m.ref"
    mdt3d = create_mt3d(ftl_path= r"K:\mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_Ash0\mf_aw_ex.ftl",
                        mt3d_name='test',
                        mt3d_ws=r"C:\Users\MattH\Downloads\test_mt3d",
                        ssm_crch=flopy.utils.Util2d.load_txt((smt.rows, smt.cols), rch_path, float, '(FREE)'),
                        ssm_stress_period_data={0:get_ssm_stress_period_data()},
                        sft_spd={0:get_sft_stress_period_data()},
                        **get_default_mt3d_kwargs())
    mdt3d.write_input()
    mdt3d.run_model()
