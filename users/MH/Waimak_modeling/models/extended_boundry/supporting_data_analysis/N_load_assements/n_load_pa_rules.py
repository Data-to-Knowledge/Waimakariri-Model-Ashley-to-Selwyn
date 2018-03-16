# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 15/03/2018 11:59 AM
"""

from __future__ import division
from core import env
import numpy as np
import pandas as pd
import geopandas as gpd
from core.spatial.vector import spatial_overlays
import os
from shapely.errors import TopologicalError

# look at dissolve in geopandas as a possible aggregation.
# farm id could be a good aggregation device

def simplify_orange_red(final_load):
    # simplify the data
    final_load = final_load.copy().dissolve(by='FARM_ID_DO', aggfunc='first')  # assumes only one landuse per farm (which is true for most (~7000 out of ~8000 shapes)
    mapper = final_load.loc[:, 'luscen_cat'].to_dict()
    keep_vars = ['geometry', 'FARM_ID_DOC', 'Farm_size_Ha',
                 'Irrig_Ha', 'Steep_Ha', 'Pstr_Ha', 'Crop_Ha', 'WF_Ha',
                 'ExtraWF_Ha', 'ExtraIrr_Ha', 'Overlap_area', 'New_NLoss_PC5', 'mean_nloss',
                 'BAU_PC5_diff_woIrrig']
    orange_red = gpd.GeoDataFrame.from_file(
        r"P:\Groundwater\Waimakariri\Landuse\Shp\Results_New_WF_rule_May17.gdb\Results_New_WF_rule_May17.gdb",
        layer='BAU_PC5_diff_woIrrig_180315')
    orange_red = orange_red.dissolve(by='FARM_ID_DOC', aggfunc=np.max)
    orange_red = orange_red.drop([' ', 'DOC', 'MASKED'])
    orange_red = orange_red.reset_index()
    orange_red.loc[:,'landuse'] = orange_red.loc[:, 'FARM_ID_DOC']
    orange_red = orange_red.replace({'landuse': mapper})
    orange_red.loc[~np.in1d(orange_red.landuse, list(set(mapper.values()))), 'landuse'] = 'Lifestyle'

    return orange_red

def _count_unique(x):
    return len(set(x))

def create_farm_scale_data(catchments, outdir):
    """

    :param catchments: dictionary of name, shapefile path
    :param red_orange_path:
    :param final_load_path:
    :param outdir:
    :return:
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    raw_out_dir = os.path.join(outdir, 'raw_data')
    if not os.path.exists(raw_out_dir):
        os.makedirs(raw_out_dir)
    final_load = gpd.GeoDataFrame.from_file(r"P:\Groundwater\Waimakariri\Landuse\N results Current Pathways.gdb",
                                            layer='Final_loadsCMPGMP140218')
    orange_red = simplify_orange_red(final_load)
    pro_rata_keys = ['pa_wforage_ha', 'pa_irrigation_ha', 'irrigated_area_ha',
                     'non_irrigated_area_ha', 'new_PA_irrigation_load', 'new_PA_forage_load',
                     'new_PA_total_load']
    outdata = pd.DataFrame(index=catchments.keys())
    for name, catchment in catchments.items():
        print(name)
        try:
            # gmp, cmp, and land_use
            temp_load = spatial_overlays(final_load, catchment)
            temp_load.loc[:, 'area_in_catch'] = temp_load.area/10000
            temp_load.loc[:, 'cmp_nload_total'] = temp_load.nload_cmp * temp_load.area_in_catch
            outdata.loc[name, 'cmp_nload_kg'] = temp_load.cmp_nload_total.sum()
            temp_load.loc[:, 'gmp_nload_total'] = temp_load.nload_gmp * temp_load.area_in_catch
            outdata.loc[name, 'gmp_nload_kg'] = temp_load.gmp_nload_total.sum()
            landuse_types = {u'',
                             u'Arable',
                             u'DairyFarm',
                             u'DairySupport',
                             u'Forest-Tussock',
                             u'Horticulture',
                             u'Lifestyle',
                             u'NotFarm',
                             u'Other-inclGolf',
                             u'Pigs',
                             u'Sheep-Beef-Deer',
                             u'SheepBeef-Hill',
                             u'Unknown'}
            grouped = temp_load.groupby('luscen_cat').aggregate({'area_in_catch': np.sum,
                                                                 'gmp_nload_total': np.sum,
                                                                 'cmp_nload_total': np.sum})
            for lt in landuse_types:
                try:
                    outdata.loc[name, '{}_area_ha'.format(lt)] = grouped.loc[lt, 'area_in_catch']
                    outdata.loc[name, '{}_gmp_load_kg'.format(lt)] = grouped.loc[lt, 'gmp_nload_total']
                    outdata.loc[name, '{}_cmp_load_kg'.format(lt)] = grouped.loc[lt, 'cmp_nload_total']
                except KeyError:
                    pass
            temp_data = spatial_overlays(orange_red, catchment)
            temp_data.loc[:, 'area_in_catch'] = temp_data.area/10000

            # convert to a dataframe
            temp_data = pd.DataFrame(temp_data.drop('geometry', axis=1))
            assert not temp_data.index.duplicated().any(), 'should not have duplicate entries from {}'.format(name)

            # create irrigated area
            temp_data.loc[:, 'irrigated_area_ha'] = temp_data.loc[:, 'Irrig_Ha']

            # create non-irrigated area
            temp_data.loc[:, 'non_irrigated_area_ha'] = temp_data.Farm_size_Ha - temp_data.irrigated_area_ha

            # create pa irrigation (sum of current pa irrigation and the potential)
            temp_data.loc[:, 'pa_irrigation_ha'] = temp_data.loc[:, 'ExtraIrr_Ha']

            # create new PA irrigation N load
            temp_data.loc[:, 'new_PA_irrigation_load'] = (temp_data.loc[:, 'New_NLoss_PC5'] - temp_data.loc[:,
                                                                                                 'BAU_PC5_diff_woIrrig']) * temp_data.Farm_size_Ha

            # create pa winter forrage (sum of current forrage and the potential)
            temp_data.loc[:, 'pa_wforage_ha'] = temp_data.loc[:, 'ExtraWF_Ha']

            # create PA forrage N load
            temp_data.loc[:, 'new_PA_forage_load'] = (temp_data.loc[:, 'BAU_PC5_diff_woIrrig'] - temp_data.loc[:,
                                                                                                    'mean_nloss']) * temp_data.Farm_size_Ha

            # create total PA N load
            temp_data.loc[:, 'new_PA_total_load'] = temp_data.loc[:, 'new_PA_forage_load'] + temp_data.loc[:,
                                                                                                     'new_PA_irrigation_load']

            # create pro-rataed data for above
            for key in pro_rata_keys:
                temp_data.loc[:, '{}_prorata'.format(key)] = temp_data.loc[:, key] * (temp_data.loc[:, 'area_in_catch'] /
                                                                                      temp_data.loc[:, 'Farm_size_Ha'])
            size_key = 'Farm_size_Ha'
            temp_data.loc[:, 'property_size_class'] = 'none'
            temp_data.loc[(temp_data.loc[:,size_key] > 0) & (temp_data.loc[:,size_key] <= 10), 'property_size_class'] = '0-10 ha'
            temp_data.loc[(temp_data.loc[:,size_key] > 10) & (temp_data.loc[:,size_key] <= 20), 'property_size_class'] = '10-20 ha'
            temp_data.loc[(temp_data.loc[:,size_key] > 20) & (temp_data.loc[:,size_key] <= 50), 'property_size_class'] = '20-50 ha'
            temp_data.loc[(temp_data.loc[:,size_key] > 50), 'property_size_class'] = '>50 ha'

            # save raw data
            temp_data.to_csv(os.path.join(raw_out_dir, name + '.csv'))

            # do some groupbys and add summary data (by property size – e.g, 0-10, 10-20, 20-50, >50 ha properties) and save data
            grouped = temp_data.groupby('landuse').aggregate(np.sum)
            grouped.loc['sum', :] = grouped.sum(axis=0)  # todo check
            grouped.to_csv(os.path.join(outdir, name + '_grouped_landuse.csv'))
            grouped = temp_data.groupby('property_size_class').aggregate(np.sum)
            grouped.loc['sum', :] = grouped.sum(axis=0)  # todo check
            grouped.to_csv(os.path.join(outdir, name + '_grouped_size.csv'))
        except TopologicalError as val:
            print(val)
    outdata.to_csv(os.path.join(outdir, 'load_overviews.csv'))


if __name__ == '__main__':
    cments = {}
    cments['whole_estuary'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\cam_bramleys_s.shp")
    cments['ashely_gorge'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\ashley_sh1.shp")
    cments['ashely_sh1'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\ash_est_all.shp")
    cments['saltwater_end'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\saltwater_end_s.shp")
    cments['waikuku_end'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\waikuku_end_s.shp")
    cments['taranaki_end'] = gpd.read_file(r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\capture_zones_particle_tracking\source_zone_polygons\likely\taranaki_end_s.shp")
    cments['loburn_fan'] = gpd.GeoDataFrame.from_file(r'\\gisdata\ProjectArchive\SCI\2015_2016\EMG\2015_2016\Waimakariri\Groundwater\BAU scenario assessment\N load analysis.gdb', layer='LoburnFanCatchment')
    #TODO SOMETHING WRONG WITH SWAZES tried vclean didn't work
    swazes = gpd.read_file(r"P:\Groundwater\Waimakariri\Surface water\Shp\Proposed SWAZ boundaries 120218\Proposed_SWAZ_Catchments.shp")
    #swazes = swazes.dissolve(by='RiverName', aggfunc='first').reset_index()
    for nm in swazes.loc[:,'RiverName']:
        cments[nm] = swazes.loc[swazes.RiverName==nm].reset_index()
    create_farm_scale_data(catchments = cments, outdir=r"C:\Users\MattH\Downloads\test_nload_stuffs")