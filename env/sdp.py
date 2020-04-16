"""
 Author: Matt Hanson
 Created: 3/11/2020 11:05 AM
 """

from __future__ import division
import os
import warnings

# todo document this when it is all said and done and set to None... and document in scripts section.
sdp_required = r"D:\Waimakariri_model_input_data\required"

assert sdp_required is not None, 'you need to set sdp_required, this is the path to the waimak data catalog'
assert os.path.exists(sdp_required), 'sdp_required does not exist, it must exist'

sdp_recommended = r"D:\Waimakariri_model_input_data\recommended"
assert sdp_recommended is not None, 'you need to set sdp_recommended, this is the path to the waimak data catalog'
assert os.path.exists(sdp_recommended), 'sdp_recommended does not exist, it must exist'

temp_file_dir =r"C:\Users\Matt Hanson\Downloads\temp_waimak_files"
assert temp_file_dir is not None, 'temp file dir must be set'
assert os.path.exists(temp_file_dir), 'temp file dir must exist'

loaded_model_realisation_dir = r"C:\Users\Matt Hanson\Downloads\temp_loaded_realisations"
assert loaded_model_realisation_dir is not None, 'loaded model realisation dir must be set'
assert os.path.exists(loaded_model_realisation_dir), 'loaded model realisation dir must exist'

pc5_ftl_repo = None # note the pc5 ftl, pc5 + eyre mar, base ftl repos are compressed in the recommended folder.
if pc5_ftl_repo is None:
    warnings.warn('pc5_ftl_repo is not set  in env.sdp.  This could cause problems for some transport modelling')

log_dir = None

if log_dir is None:
    warnings.warn('log directory is not set in sdp, this may cause exceptions in some model runs')

if not os.path.exists(temp_file_dir):
    os.makedirs(temp_file_dir)

if not os.path.exists(loaded_model_realisation_dir):
    os.makedirs(loaded_model_realisation_dir)