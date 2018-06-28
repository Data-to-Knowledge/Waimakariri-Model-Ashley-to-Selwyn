"""
Author: matth
Date Created: 28/06/2018 11:13 AM
"""

from __future__ import division
from core import env
import flopy_mh as flopy
from glob import glob
import numpy as np
import os
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt


def export_compaired_drawdown(layer, scen, outdir, percentiles):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    percentiles = np.atleast_1d(percentiles)

    paths = glob(r"D:\mh_waimak_models\stocastic_forward\*\*{}\*.hds".format(scen))
    dds = []
    for path in paths:
        try:
            base = flopy.utils.HeadFile(path.replace(scen,'current')).get_data((0,0),mflay=layer)
            base[base > 1e20] = np.nan
            new = flopy.utils.HeadFile(path).get_data((0,0),mflay=layer)
            new[new > 1e20] = np.nan
            dds.append((base-new)[np.newaxis])
        except IndexError: # avoids models which did not converge
            continue
    dds = np.concatenate(dds,axis=0)
    for per in percentiles:
        out = np.nanpercentile(dds,per, axis=0)
        smt.array_to_raster(os.path.join(outdir,'percentile_{}.tif'.format(per)),out,layer)

if __name__ == '__main__':
    export_compaired_drawdown(layer=0,scen='full_abs',outdir=r"T:\Temp\temp_gw_files\test_dds",percentiles=50)