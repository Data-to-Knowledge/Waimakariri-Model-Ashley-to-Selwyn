# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 7/09/2017 2:13 PM
"""

from __future__ import division
from model_setup import *
from data_extraction import *
from model_bc_data import *
from waimak_extended_boundry.model_run_tools.metadata_managment.convergance_check import *
from waimak_extended_boundry.model_run_tools.metadata_managment.cwms_index import *
from waimak_extended_boundry.model_run_tools.metadata_managment.transfer_readme import *
from data_extraction.data_at_wells import \
    _get_kstkpers

from data_extraction.data_from_streams import \
    _get_sw_samp_pts_dict, _get_flux_flow_arrays

# todo manage these import!!! they're hella tricky