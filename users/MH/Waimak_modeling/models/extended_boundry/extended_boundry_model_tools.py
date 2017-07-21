"""
Author: matth
Date Created: 22/06/2017 5:05 PM
"""

from __future__ import division
from core import env
from users.MH.Model_Tools.ModelTools import ModelTools
from users.MH.Waimak_modeling.supporting_data_path import sdp, temp_file_dir
from osgeo.gdal import Open as gdalOpen
import numpy as np
from users.MH.Waimak_modeling.model_tools.get_str_rch_values import get_ibound as gib
import pandas as pd
import time
import os

layers, rows, cols = 11, 364, 365

_mt = ModelTools('ex_bd_va', sdp='{}/ex_bd_va_sdp'.format(sdp), ulx=1512162.53275, uly=5215083.5772, layers=layers,
                 rows=rows, cols=cols, grid_space=200, temp_file_dir=temp_file_dir, base_mod_path=None)


def _elvdb_calc():
    top = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/tops.tif".format(sdp)).ReadAsArray()
    top[np.isclose(top, -3.40282306074e+038)] = 0

    # quickly smooth the dem near the top of the waimakariri (that gorge)
    top_smooth = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/layering/DEM_smoother.tif".format(sdp)).ReadAsArray()
    top_smooth[np.isclose(top_smooth, -3.40282306074e+038)] = 0
    top += top_smooth/2 #todo this is a first pass could be better


    # bottoms
    # layer 0
    bot0 = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/layering/base_layer_1.tif".format(sdp)).ReadAsArray()
    idx = np.isclose(bot0, -3.40282306074e+038)
    bot0[idx] = top[idx] - 10 #set nan values to 10 m thick.  all these are out of the no-flow bound

    num_layers = 9
    layer_names =[
        'ricc',
        'brom',
        'lin',
        'heath',
        'bur',
        'shirl',
        'wainoni',
        'deep_divide1',
        'deep_divide2'
    ]
    thicknesses = {
        # these values came from the median thickness from the GNS model in the confined zone
        # (defined by the extent of the ashely and chch formation)  this is a bit loose but decent given time
        'ricc': 15,
        'brom': 15,
        'lin': 25,
        'heath': 15,
        'bur': 15,
        'shirl': 15,
        'wainoni': 20,
        # seems to be target (though with minimal data) in two groups 1 from 20-30 the other from 40-80
        'deep_divide1': 40, #this one was carefully thought about to ensure that there was the potential to grab targets (assuming no changes in the above layers)
        'deep_divide2': 40 #todo this captures our deepest wells but means no targets in teh bottom of the model.

    }
    outdata = np.zeros((len(layer_names)+3,_mt.rows,_mt.cols))
    outdata[0] = top
    outdata[1] = bot0

    for i,name in enumerate(layer_names):
        layer = i+2
        thickness = thicknesses[name]
        outdata[layer] = outdata[layer-1] - thickness

    # bottom of the model
    basement = _get_basement()
    thickness = outdata[-2] - basement
    thickness[thickness < 20] = 20
    outdata[-1] = outdata[-2] - thickness

    waihora = _mt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/te_waihora.shp".format(_mt.sdp), 'ID', True)
    idx = np.isfinite(waihora)
    outdata[0][idx] = 1.5

    temp_path = env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and optimisation/stream_drn_values/ashley estuary.shp")
    ash_est = smt.shape_file_to_model_array(temp_path, 'ID', alltouched=True)
    outdata[0][np.isfinite(ash_est)] = 1.2

    # hackish fix to keep stream at or below ground surface
    outdata[0,57,225] +=1
    return outdata

def _get_basement():
    basement = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/basement2.tif".format(sdp)).ReadAsArray()
    basement[np.isclose(basement, 9999999)] = np.nan
    basement = basement[1:,:]
    basement2 = np.concatenate((basement[:,:], np.repeat(basement[:, -1][:, np.newaxis], 33, axis=1)), axis=1)
    return basement2


def _no_flow_calc():
    no_flow = np.zeros((_mt.rows,_mt.cols))
    outline = _mt.shape_file_to_model_array("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/new_active_domain.shp".format(sdp),'DN',True)
    no_flow[np.isfinite(outline)] = 1

    no_flow = np.repeat(no_flow[np.newaxis,:,:],_mt.layers, axis=0)
    # convert shapefile to array and set all to 1 then add to the active boundry (set all <0 to 1)


    # propogate constant head of the sea down as active cells (handle mis-match of outline and sea values
    sea = _mt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/coastal_constant_heads.shp".format(_mt.sdp),'ID',True)
    for i in range(1,_mt.layers):
        no_flow[i][np.isfinite(sea)] = 1

    # set no flow from boundry
    tops = _elvdb_calc()[0:_mt.layers]
    if len(tops) != _mt.layers:
        raise ValueError('wrong number of layers returned')
    basement = np.repeat(_get_basement()[np.newaxis,:,:], _mt.layers, axis=0)
    buffer = np.array([0,5,5,5,5,5,5,5,10,10,10]) # set all within buffer to no flow
    no_flow[basement+buffer[:, np.newaxis, np.newaxis] >= tops] = 0

    # fix the pockets of no-flow inside the model, and general domain
    no_flow_clip = np.zeros((smt.rows,smt.cols))*np.nan
    for i in range(smt.layers):
        shp_path = (env.sci("Groundwater/Waimakariri/Groundwater/Numerical GW model/Model build and "
                                      "optimisation/stream_drn_values/no flow rasters/noflowpockets_"
                                      "layer{:02d}.shp".format(i)))
        if os.path.exists(shp_path):
            temp = smt.shape_file_to_model_array(shp_path,'Id',True)
            no_flow_clip[np.isfinite(temp)] = 0

        no_flow[i][np.isfinite(no_flow_clip)] = 0


    #set constant heads  to -1
    constant_heads = _get_constant_heads()
    no_flow[np.isfinite(constant_heads)] = -1
    no_flow[(np.isclose(constant_heads, -999))] = 1

    return no_flow


def _get_constant_heads():

    outdata = np.zeros((_mt.layers, _mt.rows, _mt.cols))*np.nan

    # sea surface (north/south wai
    first_sea_val = 0.5
    rest_sea_val = 0
    sea = _mt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/coastal_constant_heads.shp".format(_mt.sdp),'ID',True)
    sea[np.isfinite(sea)] = rest_sea_val

    # set east most value to first sea for the north part of the model
    for i in range(0,300):
        if all(np.isnan(sea[i])):
            continue
        j = np.where(np.isfinite(sea[i]))[0].min()
        sea[i,j] = first_sea_val

    temp = sea[300:,:]
    for j in range(309):
        if all(np.isnan(temp[:,j])):
            continue

        i = np.where(np.isfinite(temp[:,j]))[0].min()
        sea[i+300,j] = first_sea_val

    outdata[0][np.isfinite(sea)] = sea[np.isfinite(sea)]
    # propogate constant head down as active cells
    for i in range(1,_mt.layers):
        outdata[i][np.isfinite(sea)] = -999

    tops = _elvdb_calc()[0:_mt.layers]
    if len(tops) != _mt.layers:
        raise ValueError('wrong number of layers returned')
    basement = np.repeat(_get_basement()[np.newaxis,:,:], _mt.layers, axis=0)
    buffer = np.array([5,5,5,5,5,5,5,5,10,10,10]) # set all within buffer to no flow
    outdata[basement+buffer[:, np.newaxis, np.newaxis] >= tops] = np.nan

    return outdata


smt = ModelTools(
    'ex_bd_va', sdp='{}/ex_bd_va_sdp'.format(sdp), ulx=1512162.53275, uly=5215083.5772, layers=layers, rows=rows,
    cols=cols, grid_space=200, no_flow_calc=_no_flow_calc, temp_file_dir=temp_file_dir, elv_calculator=_elvdb_calc,
    base_mod_path=None
)

# quick versioning
model_version = 'a'

if model_version == 'a':
    smt.model_version = 'a'
    smt.sfr_version, smt.seg_v, smt.reach_v = 1, 1, 1
    smt.k_version = 1
    smt.wel_version = 1



if __name__ == '__main__':
   start = time.time()

   print('took {} seconds'.format((time.time()-start)))
