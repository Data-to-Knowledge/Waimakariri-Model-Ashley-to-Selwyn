"""
Author: matth
Date Created: 7/09/2017 4:19 PM
"""

from __future__ import division
import os
from waimak_extended_boundry import smt
from waimak_extended_boundry.model_run_tools.model_setup.base_modflow_wrapper import mod_gns_model
from waimak_extended_boundry.model_run_tools.model_bc_data.wells import  get_race_data, get_full_consent
from waimak_extended_boundry.model_run_tools.data_extraction.data_at_wells import hds_no_data, get_hds_file_path
import flopy_mh as flopy
import numpy as np
from waimak_extended_boundry.model_run_tools.model_bc_data.drn_data import get_drn_no_ncarpet_spd
from stream_depletion_sdp import starting_heads_dir


def get_sd_starting_hds(model_id, sd_version):
    if sd_version == 'sd7':
        hds = get_starting_heads_sd7(model_id)
    elif sd_version == 'sd30':
        hds = get_starting_heads_sd30(model_id)
    elif sd_version == 'sd150':
        hds = get_starting_heads_sd150(model_id)
    elif sd_version == 'grid':
        hds = _get_no_pumping_ss_hds(model_id)
    else:
        raise ValueError('unexpected argument for sd_version: {}'.format(sd_version))

    return hds

def get_sd_well_list(model_id):
    cav = get_full_consent(model_id, missing_sd_wells=True) # note this is identical to max rate well list
    cav = cav.loc[(cav.type=='well') & (cav.zone == 'n_wai') & (cav.flux < 0)]
    well_list = list(set(cav.index))
    return well_list

def get_starting_heads_sd150(model_id):
    hds = _get_no_pumping_ss_hds(model_id)
    return hds


def get_starting_heads_sd30(model_id):
    hds = _get_no_pumping_ss_hds(model_id)
    return hds


def get_starting_heads_sd7(model_id):
    hds = _get_no_pumping_ss_hds(model_id)
    return hds


def _get_no_pumping_ss_hds(model_id, recalc=False):

    pickle_path = "{}/model_{}_sd_starting_hds.txt".format(starting_heads_dir, model_id)
    if (os.path.exists(pickle_path)) and (not recalc):
        hds = np.loadtxt(pickle_path)
        if hds.shape != (smt.layers, smt.rows, smt.cols):
            raise ValueError('unexpected shape for hds {}, expected {}'.format(hds.shape, (smt.layers, smt.rows, smt.cols)))
        return hds
    dirpath = "{}/forward_supporting_models/base_str_dep".format(smt.sdp)  # model Id is added in import gns model
    drns = get_drn_no_ncarpet_spd(model_id)
    well = {0: smt.convert_well_data_to_stresspd(get_race_data(model_id))}
    m = mod_gns_model(model_id, 'base_heads_for_str_dep', dir_path=dirpath, safe_mode=False, well=well, drain={0:drns})
    m.write_name_file()
    m.write_input()
    suc, buff = m.run_model()
    if not suc:

        raise ValueError('starting heads model for model_id: {} did not run'.format(model_id))
    hds_path = get_hds_file_path(m=m)
    hds = flopy.utils.HeadFile(hds_path).get_alldata(nodata=hds_no_data)[-1, :, :, :]
    if hds.shape != (smt.layers, smt.rows, smt.cols):
        raise ValueError('unexpected shape for hds {}, expected {}'.format(hds.shape, (smt.layers, smt.rows, smt.cols)))
    np.savetxt(pickle_path, hds)
    return hds


def get_ss_sy(ss_sy_version=1, return_description=False):
    """
    a way of keeping track of the versioning
    :param ss_sy_version: numeric
    :param return_description: bool if True do not return the ss_sy, but instead return the description for the ss_sy
                               version for writing into a readme file
    :return: ss, sy or description
    """
    if ss_sy_version == 1:
        sy = np.zeros(smt.model_array_shape) + 0.1
        ss = np.zeros(smt.model_array_shape)
        ss[0, :, :] = 1.6E-3
        ss[1:, :, :] = 5.6E-5
        description = ('the first pass storage, created from average storage conditions, but not considering '
                      'the implication of flooded cells\n sy: 0.1 for all layers\n'
                      'ss: 1.6E-3 for layer 1 and 5.6 E-5 for all other layers')
    elif ss_sy_version == 2:
        sy = smt.get_empty_model_grid() + 0.01
        ss = smt.get_empty_model_grid() + 1.0E-6
        description = 'the lower of the sensitivity anlaysis\n ss: 1e-6, sy: 0.01'
    elif ss_sy_version == 3:
        sy = smt.get_empty_model_grid() + 0.05
        ss = smt.get_empty_model_grid() + 1.0E-5
        description = 'the median values for the sensitivity analysis\n ss: 1e-5, sy: 0.05'
    elif ss_sy_version == 4:
        sy = smt.get_empty_model_grid() + 0.1
        ss = smt.get_empty_model_grid() + 5E-4
        description = 'the high values for the sensitivity analysis\n ss: 5e-4, sy: 0.1'
    else:
        raise ValueError('unexpected ss_sy_version')

    if return_description:
        return description
    else:
        return ss, sy


if __name__ == '__main__':
    test = get_sd_well_list('NsmcBase')
    print(len(test))