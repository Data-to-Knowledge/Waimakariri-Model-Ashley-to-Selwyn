# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 9/05/2018 8:37 AM
"""

from __future__ import division
import geopandas as gpd

if __name__ == '__main__':
    final_load = gpd.GeoDataFrame.from_file(r"P:\Groundwater\Waimakariri\Landuse\N results Current Pathways.gdb",
                                            layer='Final_loadsCMPGMP140218')
    final_load.groupby('luscen_cat')['nload_gmp'].describe()