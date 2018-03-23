# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 23/03/2018 2:27 PM
"""

from __future__ import division
from core import env
from multiple_source_delineation import create_amalgimated_source_protection_zones, split_netcdfs
import os


if __name__ == '__main__':
    model_id = 'NsmcReal{:06d}'.format()
    recalc=True
    recalc_backward_tracking=False
    base_outdir = r"D:\mh_waimak_models\private_domestic_supply"
    print('running for AshOpt')
    outdir = os.path.join(base_outdir, model_id)
    create_amalgimated_source_protection_zones(model_ids=['AshOpt'], run_name='AshOpt_private_wells',
                                               outdir=outdir, num_it=1,
                                               recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(outdir)
