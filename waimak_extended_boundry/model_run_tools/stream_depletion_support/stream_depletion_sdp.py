"""
 Author: Matt Hanson
 Created: 16/04/2020 11:51 AM
 """
import warnings

#todo look through documentation

base_sd_dir = None

if base_sd_dir is None:
    warnings.warn('base sd is not set in '
                  'waimak_extended_boundry/model_run_tools/stream_depletion_support/stream_depletion_sdp.py, '
                  'this will cause problems with stream depletion modelling'
                  )

starting_heads_dir = None

if starting_heads_dir is None:
    warnings.warn('starting head dir is not set in '
                  'waimak_extended_boundry/model_run_tools/stream_depletion_support/stream_depletion_sdp.py, '
                  'this will cause problems with stream depletion modelling'
                  )