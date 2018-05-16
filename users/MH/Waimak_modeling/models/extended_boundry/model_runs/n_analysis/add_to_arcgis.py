# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 17/04/2018 2:10 PM
"""

if __name__ == '__main__':
    import arcpy

arcpy.env.workspace = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\n_red_mt3d"
import os

workspace = arcpy.env.workspace
ncList = arcpy.ListFiles("*.tif")
symbologyLayer = "Symbology.lyr" # symbology layer is in same file as teh arcpy workspace above

for i, nc in enumerate(ncList):
    if 'reason' in nc:
        symbol = 'reason_template.lyr'
    else:
        symbol = 'reduction_template.lyr'
    arcpy.MakeRasterLayer_management(os.path.join(workspace,nc), nc.replace('.tif',''))
    arcpy.ApplySymbologyFromLayer_management(nc.replace('.tif',''), os.path.join(workspace,symbol))

# to add shapefiles this does not work
arcpy.env.workspace = r"P:\Groundwater\Waimakariri\Groundwater\Numerical GW model\Model simulations and results\ex_bd_va\n_results\interzone_n_results\receptors_with_data"
import os

workspace = arcpy.env.workspace
ncList = arcpy.ListFiles("*.shp")
symbologyLayer = 'layer_with_table.lyr' # symbology layer is in same file as teh arcpy workspace above

for i, nc in enumerate(ncList):
    arcpy.MakeFeatureLayer_management(os.path.join(workspace,nc), nc.replace('.shp',''))


mxd = arcpy.mapping.MapDocument("CURRENT")
layers = arcpy.mapping.ListLayers(mxd)
df = arcpy.mapping.ListDataFrames(mxd)[0]
for lay in layers:
    if lay.name == 'layer for using':
        sourceLayer = lay
    else:
        continue
for lay in layers:
    if lay.name == 'layer for using':
        continue
    arcpy.mapping.UpdateLayer(df, lay, sourceLayer, True)
