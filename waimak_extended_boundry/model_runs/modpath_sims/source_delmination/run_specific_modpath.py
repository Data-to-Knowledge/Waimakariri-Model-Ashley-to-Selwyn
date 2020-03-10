# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 23/03/2018 2:27 PM
"""

from __future__ import division
from single_zone_delination import create_single_zone_indexs,create_zones, split_netcdfs
import os


if __name__ == '__main__':
    model_id = 'NsmcReal{:06d}'.format(491)
    recalc=False
    recalc_backward_tracking=False
    indexes = create_single_zone_indexs()
    base_outdir = r"D:\mh_waimak_models\single_source_zones"
    print('running for {}'.format(model_id))
    outdir = os.path.join(base_outdir, model_id)
    create_zones(model_ids=[model_id], run_name='{}_single_sources'.format(model_id),
                 outdir=outdir, root_num_part=3,
                 indexes=indexes, recalc=recalc, recalc_backward_tracking=recalc_backward_tracking)
    split_netcdfs(outdir)

