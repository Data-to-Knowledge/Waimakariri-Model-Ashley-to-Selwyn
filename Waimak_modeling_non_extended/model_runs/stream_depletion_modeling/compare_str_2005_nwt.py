"""
Author: matth
Date Created: 29/05/2017 3:04 PM
"""

from __future__ import division

from users.MH.Waimak_modeling.model_tools.get_str_flows import streamflow_for_kskps
from users.MH.Waimak_modeling.supporting_data_path import sdp, base_mod_dir

m2005_path = "{}/base_model_runs/base_ss_mf/base_SS".format(sdp)

mnwt_path = '{}/s_explore_transient/perc_50_ss_0.000169836616135_sy_0.000472232863272/perc_50_ss_0.000169836616135_sy_0.000472232863272'.format(base_mod_dir)

m2005 = streamflow_for_kskps(m2005_path, (0,0))
mnwt = streamflow_for_kskps(mnwt_path, (0, 0))

print ((m2005-mnwt)/m2005*100)