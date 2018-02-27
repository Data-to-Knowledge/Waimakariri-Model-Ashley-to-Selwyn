# -*- coding: utf-8 -*-
"""
Author: mattH
Date Created: 30/08/2017 2:04 PM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import os
import geopandas as gpd
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.NSMC_inputs.recharge_index_array import get_rch_index_array
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.NSMC_inputs.gen_well_csv_for_NSMC import get_well_NSMC_base, get_template_data
from users.MH.Waimak_modeling.models.extended_boundry.extended_boundry_model_tools import smt
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.NSMC_inputs.parameter_file import create_parameter_file, get_param_table
from shutil import copyfile
from users.MH.Waimak_modeling.models.extended_boundry.supporting_data_analysis.NSMC_inputs.NSMC_weights import get_NSMC_weights
def data_to_cath(base_dir):
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    # rch_array
    rch_idx = get_rch_index_array()
    np.savetxt('{}/rch_idx.txt'.format(base_dir),rch_idx)

    #  rch_points white delim
    points = gpd.read_file("{}/m_ex_bd_inputs/shp/rch_pilot_points_option2.shp".format(smt.sdp))
    # groups: 1: dryland, 2: confined, 3: selwyn_irr, 4: waimak_irr
    points.loc[:,'name'] = ['rch_ppt_{:02d}'.format(e) for e in range(len(points.index))]
    points.loc[:, 'x'] = [e.x for e in points.geometry]
    points.loc[:, 'y'] = [e.y for e in points.geometry]
    points.loc[:, 'value'] = ['${}$'.format(e) for e in points.loc[:,'name']]
    points2 = points.loc[:,['name','x','y','group']]
    points2.to_csv('{}/rch_ppts.txt'.format(base_dir),header=False, sep=' ',index=False,mode='a')
    out_file = '{}/rch_ppts.tpl'.format(base_dir)
    points3 = points.loc[:,['name','x','y','group','value']]
    with open(out_file,'w') as f:
        f.write('ptf $\n')
    points3.to_csv(out_file,header=False, sep=' ',index=False,mode='a')

    # well template
    well_temp = get_template_data()
    outpath = '{}/wel_modifiers.tpl'.format(base_dir)
    with open(outpath,'w') as f:
        f.write('ptf $\n')
    well_temp.to_csv(outpath,header=False, sep=' ',mode='a')

    # well base csv
    all_wells =  get_well_NSMC_base()
    all_wells.to_csv('{}/wel_pkg_base_data.csv'.format(base_dir))

    # parameter file this is not correct yet cath is finishing it
    create_parameter_file('{}/nsmc_par.unc'.format(base_dir))

    # script?
    copyfile(r'C:\Users\MattH\OneDrive - Environment Canterbury\Ecan_code4\Ecan.Science.Python.Base\users\MH\Waimak_modeling\models\extended_boundry\supporting_data_analysis\NSMC_inputs\well_pkg_adjust.py','{}/well_pkg_adjust.py'.format(base_dir))

    # parameter table:
    parm = get_param_table()
    parm.to_csv('{}/parameter_table.csv'.format(base_dir))

    #NSMC_weight
    nsmc_weights = get_NSMC_weights()
    nsmc_weights.to_csv('{}/nsmc_weights.txt'.format(base_dir), header=False, sep=' ')


if __name__ == '__main__':
    data_to_cath(r"C:\Users\MattH\Desktop\NSMC_data_to_cath_check")