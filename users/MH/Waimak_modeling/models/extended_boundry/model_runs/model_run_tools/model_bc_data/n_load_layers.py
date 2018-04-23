# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 18/04/2018 11:20 AM
"""

from __future__ import division
from core import env
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
import numpy as np
import os
import flopy


def get_gmp_con_layer(recalc=False):
    pickle_path = "{}/gmp_n_conc.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_new_cmp_con(recalc=False):
    # should be similar to the get original cmp layer, but I am unsure how the re-sampling was handled for that layer
    pickle_path = "{}/cmp_n_conc_resampled.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nconc_cmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_pc5pa_additonal_load(recalc=False):
    pickle_path = "{}/pc5pa_additional_n_load.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb"
    layer = 'BAU_PC5_diff_woIrrig_180315'

    outdata = smt.geodb_to_model_array(n_load_path, layer, attribute='New_NLoss_PC5', alltouched=True,
                                       area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_gmp_load_raster(recalc=False):
    pickle_path = "{}/gmp_n_load.txt".format(smt.pickle_dir)
    if (os.path.exists(pickle_path)) and (not recalc):
        outdata = np.loadtxt(pickle_path)
        return outdata

    n_load_path = env.sci('Groundwater\\Waimakariri\\Groundwater\\Numerical GW model\\Model simulations and '
                          'results\\Nitrate\\NloadLayers\\CMP_GMP_PointSources290118_nclass.shp')
    outdata = smt.shape_file_to_model_array(n_load_path, attribute='nload_gmp', alltouched=True,
                                            area_statistics=True, fine_spacing=10, resample_method='average')
    outdata[np.isnan(outdata)] = 0
    np.savetxt(pickle_path, outdata)
    return outdata


def get_orginal_cmp_layer():
    # the layer that all of brioch's runs were done with
    rch_path = env.gw_met_data(r"mh_modeling\data_from_gns\AshOpt_medianN\AWT20180103_Ash0\AWT20180103_"
                               r"Ash0\nconc_cmp_200m.ref")

    return flopy.utils.Util2d.load_txt((smt.rows, smt.cols), rch_path, float, '(FREE)')


if __name__ == '__main__':
    if True:
        gmp = get_gmp_con_layer()
        cmp = get_new_cmp_con()
        temp = gmp / cmp
        smt.array_to_raster(
            r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\interzone_n_results\gmpcon_over_cmpcon.tif",
            temp, no_flow_layer=0)
