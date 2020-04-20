"""
 Author: Matt Hanson
 Created: 3/11/2020 11:05 AM

 The purpose of this script is to locate the larger required and recommended files associated with the model.
 These cannot be stored in the git repository as they are large (c. 264 GB).  several directories need to be set
 specifically:
    sdp_required: part of the waimakariri Model data catalog, this is essential and contains the data needed to create
                  model realisations including all of the boundary condition data this is (c. 17 GB).
                  this is frequently accessed so it is best to move it onto a faster drive.
    sdp_recommended: important data, but only needed for some of the more complex modelling or to reduce run times.
                     however this data set is large (c. 240 GB) for simpiler modelling it may be possible to set this
                     value to a blank folder, which will minimize required space.


several directories need to be set and made, but they can be empty (stuff will be written here)
    temp_file_dir: this is just a temp file dir, not much will be saved here and most will be removed programatically
                   only 1-2 GB should be sufficient for this.
    loaded_model_realisation_dir: this is a place for the model to save loaded realisations.  in this modelling the
                                  models are saved as their parameter sets.  this allows a model to be stored
                                  efficiently. The modeliling process then creates a modflow model from this and saves
                                  it to this loaded directory.  It does this because it takes c. 4-15 min to load a
                                  realisation from the parameter sets and c. 15s to load it from the saved dir. Note
                                  that saving the file to a directory takes c.  0.25 gb/model,

additionally several directories are needed for some more complex modelling:
    pc5_ftl_repo: a set of mt3d link files, this exists in recommended, but is compressed
    log_dir: a place to write log files.
 """

from __future__ import division
import os
import warnings

#todo set to None for distribution
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
