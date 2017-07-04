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

layers, rows, cols = 11, 364, 365

_mt = ModelTools('ex_bd_va', sdp='{}/ex_bd_va_sdp'.format(sdp), ulx=1512162.53275, uly=5215083.5772, layers=layers,
                 rows=rows, cols=cols, grid_space=200, temp_file_dir=temp_file_dir, base_mod_path=None)


def _elvdb_calc():  # todo, add pickle
    top = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/tops.tif".format(sdp)).ReadAsArray()
    top[np.isclose(top, -3.40282306074e+038)] = 0


    # bottoms
    # layer 0
    bot0 = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/layering/base_layer_1.tif".format(sdp)).ReadAsArray()  # todo define this
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
    thicknesses = {  #todo check verticle offset of targets, particularly vertical gradient targets
        # these values came from the median thickness from the GNS model in the confined zone
        # (defined by the extent of the ashely and chch formation)  this is a bit loose but decent given time
        'ricc': 15, #todo may need to smooth this layer inland? it may be worth making this layer non-uniform
        'brom': 15,
        'lin': 25,
        'heath': 15,
        'bur': 15, #todo this could pose a problem for the steep sections of the model layering
        'shirl': 15,
        'wainoni': 20,
        # seems to be target (though with minimal data) in two groups 1 from 20-30 the other from 40-80
        'deep_divide1': 40, #this one was carefully thought about to ensure that there was the potential to grab targets (assuming no changes in the above layers)
        'deep_divide2': 40 #todo this captures our deepest wells but means no targets in teh bottom of the model.  I could join these two... think about v gradient
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





    #todo check that no layer is less than 50% of it's neighbor #todo implement in model tools

    waihora = _mt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/te_waihora.shp".format(_mt.sdp), 'ID', True)
    idx = np.isfinite(waihora)
    outdata[0][idx] = 1.5 #todo check... set DEM to lake level

    return outdata

def _get_basement():
    basement = gdalOpen("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/basement2.tif".format(sdp)).ReadAsArray()
    basement[np.isclose(basement, 9999999)] = np.nan
    basement = basement[1:,:]
    basement2 = np.concatenate((basement[:,:], np.repeat(basement[:, -1][:, np.newaxis], 33, axis=1)), axis=1)
    return basement2


def _no_flow_calc():  # todo, add pickle if it takes long... shouldn't
    no_flow = np.zeros((_mt.rows,_mt.cols))
    # todo weird artifact here in
    outline = _mt.shape_file_to_model_array("{}/ex_bd_va_sdp/m_ex_bd_inputs/shp/new_active_domain.shp".format(sdp),'DN',True)
    no_flow[np.isfinite(outline)] = 1

    no_flow = np.repeat(no_flow[np.newaxis,:,:],_mt.layers, axis=0)
    # convert shapefile to array and set all to 1 then add to the active boundry (set all <0 to 1)

    # set no flow from boundry
    tops = _elvdb_calc()[0:_mt.layers]
    if len(tops) != _mt.layers:
        raise ValueError('wrong number of layers returned')
    basement = np.repeat(_get_basement()[np.newaxis,:,:], _mt.layers, axis=0)
    no_flow[basement >= tops] = 0

    #todo check the pockets of no-flow inside the model, and general domain

    # set constant heads  to -1 #todo implement once constant heads are defined
    #constant_heads = _get_constant_heads()
    #no_flow[np.isfinite(constant_heads)] = -1

    return no_flow



def _get_constant_heads():

    outdata = np.zeros((_mt.layers, _mt.rows, _mt.cols))*np.nan
    # sea surface (north/south wai
    sea_val = None #todo
    sea = _mt.shape_file_to_model_array(None,'ID',True) #todo build shapefile
    idx = np.isfinite(sea)
    outdata[0][idx] = sea_val

    # te waihuro
    wai_val = 1.5 #todo check
    waihora = _mt.shape_file_to_model_array("{}/m_ex_bd_inputs/shp/te_waihora.shp".format(_mt.sdp), 'ID', True)
    idx = np.isfinite(waihora)
    outdata[0][idx] = wai_val

    # todo model boundry (sw)... how many layers for this
    # return a 3d array (layer, col, row) of constant heads values and np.nan for all others.


    raise NotImplementedError


smt = ModelTools(
    'ex_bd_va', sdp='{}/ex_bd_va_sdp'.format(sdp), ulx=1512162.53275, uly=5215083.5772, layers=layers, rows=rows,
    cols=cols, grid_space=200, no_flow_calc=_no_flow_calc, temp_file_dir=temp_file_dir, elv_calculator=_elvdb_calc,
    base_mod_path=None
)



if __name__ == '__main__':
   _no_flow_calc()

   print 'done'