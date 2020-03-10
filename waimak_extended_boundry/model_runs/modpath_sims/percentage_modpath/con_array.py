# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 8/12/2017 1:28 PM
"""

from __future__ import division
import geopandas as gpd
from waimak_extended_boundry import smt
import pandas as pd
import flopy_mh as flopy
from waimak_extended_boundry import _get_wel_spd_v3

def make_con_layer(bnd_type, rch_con_array, well_data, sfr_data):
    """
    make a concentration layer for a given model realisation, this assumes no transport through surface water
    :param bnd_type: boundary type array (from particle generation)  0=rch,1=well,2==sfr, -1 is no particles
    :param rch_con_array: a (i,j) array of recharge concentrations
    :param well_data: pd.DataFrame concentration data for the injection wells with conc as the id variable
    :param sfr_data: pd.DataFrame concentration for each sfr reach with conc as the id variable
    :return: np.ndarray shape: (smt.rows,smt.cols)
    """
    outdata = smt.get_empty_model_grid()

    well_con_array = smt.df_to_array(well_data, 'conc')
    sfr_con_array = smt.df_to_array(sfr_data, 'conc')
    # populate recharge concentrations
    outdata[bnd_type==0] = rch_con_array[bnd_type==0]
    outdata[bnd_type==1] = well_con_array[bnd_type==1]
    outdata[bnd_type==2] = sfr_con_array[bnd_type==2]

    return outdata

def make_inital_sfr_dataframe(ashley, cust, eyre, waimak):
    """
    make an inital dataframe assumes that all cust and ashley tribs are identical and no stream routing
    :param ashley: concentration for the top of the ashely
    :param cust: concentration for the top of the ashely
    :param eyre: concentration for the top of the ashely
    :param waimak: concentration for the top of the ashely
    :return:
    """
    sfr = gpd.read_file(r"{}\m_ex_bd_inputs\raw_sw_samp_points\sfr\all_sfr.shp".format(smt.sdp))
    sfr = sfr.loc[:,['riv_name', 'k', 'i', 'j']]
    sfr.loc[:,'conc'] = 0
    sfr.loc[sfr['riv_name']=='waimak','conc'] = waimak
    sfr.loc[sfr['riv_name']=='eyre','conc'] = eyre
    sfr.loc[sfr['riv_name']=='cust','conc'] = cust
    sfr.loc[sfr['riv_name']=='ashley','conc'] = ashley

    return sfr


def _make_mednload_approx(bnd_type):
    """
    make an n load
    :param bnd_type:
    :return:
    """

    sfr = make_inital_sfr_dataframe(0.1,0.35,0.35,0.1)
    temp_well = pd.read_table(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\median_n_load_wells.txt",
                         delim_whitespace=True, index_col=0)
    temp_well.loc[:,['layer','row','col']] += -1
    well = _get_wel_spd_v3()
    well.loc[:,'conc'] = 0
    well.loc[temp_well.index,'conc'] = temp_well.loc[:,'conc']
    rch = flopy.utils.Util2d.load_txt((smt.rows,smt.cols),r"K:\mh_modeling\data_from_gns\run_files\nconc_cmp_200m.ref",float,'(FREE)')

    load = make_con_layer(bnd_type,rch,well,sfr)
    return load

# could try to handle stream flow routing, but this didn't work

if __name__ == '__main__':
    import os
    import numpy as np
    mp_ws = r"D:\mh_waimak_models\modpath_emulator"
    mp_name = 'NsmcBase_first_try'
    bnd_type = np.loadtxt(os.path.join(mp_ws,'{}_bnd_type.txt'.format(mp_name)))
    load = _make_mednload_approx(bnd_type)
    print'done'
