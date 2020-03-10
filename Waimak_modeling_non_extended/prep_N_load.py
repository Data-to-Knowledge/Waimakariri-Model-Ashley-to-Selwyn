"""
Author: matth
Date Created: 7/04/2017 1:19 PM
"""

from __future__ import division
import env
import numpy as np
import model_tools as mt
from copy import deepcopy
from supporting_data_path import temp_file_dir

type = 2
if type ==1:
    basedir = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model")
    outdir = '{}/Model build and optimisation/Nitrate/N load inputs/N_load_matrix'.format(basedir)
    gdb_path = '{}/Model build and optimisation/Nitrate/N load inputs/N load layers/N_Input_Layers.gdb'.format(basedir)

    recharge = mt.read_visas_matrix('{}/recharge_from_model.DAT'.format(outdir))
    recharge[recharge <= 0] = np.nan
    recharge[mt.no_flow] = np.nan
    recharge2 = deepcopy(recharge)
    recharge2[np.isnan(recharge2)] = 0
    mt.write_vistas_matrix('{}/positive_recharge.dat'.format(outdir),recharge2)

    shp_files = ['nloss_CurrentMP_adj'] #start for now , 'nloss_GMP_adj', 'nloss_PC5_orangeNAZ', 'nloss_PC5_redNAZ' ]
    attributes = ['nloss_om_adj']
    outnames = ['CMP_N_load']

    for shpfile, att, outname in zip(shp_files,attributes,outnames):
        load_data_temp = mt.geodb_to_model_array(gdb_path, shpfile, att)

        # move from kg/m2/year to g/m3 (concentration of recharge)
        load_data = (load_data_temp * 1000 / 365 /10000)/ recharge

        load_data[np.isnan(load_data)] = 0
        load_data[mt.no_flow] = 0
        load_data[load_data > 200] = 200
        fig, ax = mt.plt_matrix(load_data)
        fig.show()
        mt.write_vistas_matrix('{}/{}.dat'.format(outdir,outname), load_data)

if type == 2:
    load_data_temp = mt.shape_file_to_model_array(
        env.sci(r"Groundwater\Waimakariri\Groundwater\Groundwater Quality\Shp\N_Conc_CP_adj_Rough.shp")
        ,'N_CMP_Conc')
    load_data_temp[load_data_temp>200] = 200
    load_data_temp[np.isnan(load_data_temp)] = 0
    print(load_data_temp.shape)
    mt.write_vistas_matrix(r"{}/first_n_con_test.txt".format(temp_file_dir), load_data_temp)
