"""
Author: matth
Date Created: 1/05/2017 1:43 PM
"""

from __future__ import division
from core import env
import os
import socket

sdp = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/supporting_data_for_scripts")

results_dir = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model simulations and results")


temp_file_dir = env.transfers("Temp/temp_gw_files")
if not os.path.exists(temp_file_dir):
    os.makedirs(temp_file_dir)

_comp = socket.gethostname()

if _comp == 'HP1639':
    base_mod_dir = "C:/Users/MattH/Desktop/Waimak_modeling/python_models"
    base_mod_dir2 = base_mod_dir
elif _comp == 'DHI-Runs02':
    base_mod_dir = "D:/mattH/python_wm_runs"
    base_mod_dir2 = "E:/python_wai_models" # to handle a second drive on DHI-Runs 02
elif _comp == 'GWATER02':
    base_mod_dir = "D:/mh_waimak_models"
    base_mod_dir2 = base_mod_dir
elif _comp == 'RDSProd03':
    base_mod_dir = r"C:\Users\matth\Desktop\base_mod"
    base_mod_dir2 = base_mod_dir
elif _comp == 'Desktop':
    base_mod_dir = r"C:\matt_modelling_unbackedup"
    base_mod_dir2 = base_mod_dir
else:
    raise ValueError('unidentified machine {} no default path for models'.format(_comp))


def get_org_mod_path(model='m_strong_vert'):
    """
    a function to keep track of all our model versions
    :param model: the model to select (see keys below) defaults to m_strong_vert to ensure backwards compatability for
                  runs
    :return: path to the model (excluding extension)
    """
    paths = {
        'm_flooded': "{}/from_GNS/m_flooded/BH_OptMod_Flow".format(sdp),
        'm_strong_vert': '{}/from_GNS/m_strong_vert/native_mf/BH20170502_m2/BH20170210_bh0501'.format(sdp),
        'm_relax_vert': "{}/from_GNS/m_relax_vert/native_mf/AWFlow20170502_m3_origtxt/BH20170210_bh0501".format(sdp)
    }
    if model not in paths.keys():
        raise ValueError('model {} not defined'.format(model))

    return paths[model]
