# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# gauge_loc_v02.py
# Created on: 2016-04-07 10:09:48.00000
#   (generated by ArcGIS/ModelBuilder)
# Description:
# ---------------------------------------------------------------------------

# Set the necessary product code
# import arcinfo


# Import arcpy module
import arcpy
from arcpy import env
from arcpy.sa import *
#import ArcHydroTools as ah

# Define functions

def snap_points(points, lines, distance):

    import arcgisscripting, sys

    gp = arcgisscripting.create()

    # Load the Analysis toolbox so that the Near tool is available
    gp.toolbox = "analysis"

    # Perform the Near operation looking for the nearest line
    # (from the lines Feature Class) to each point (from the
    # points Feature Class). The third argument is the search
    # radius - blank means to search as far as is needed. The
    # fourth argument instructs the command to output the
    # X and Y co-ordinates of the nearest point found to the
    # NEAR_X and NEAR_Y fields of the points Feature Class
    gp.near(points, lines, str(distance), "LOCATION")

    # Create an update cursor for the points Feature Class
    # making sure that the NEAR_X and NEAR_Y fields are included
    # in the return data
    rows = gp.UpdateCursor(points, "", "", "NEAR_X, NEAR_Y")

    row = rows.Next()

    # For each row
    while row:
        # Get the location of the nearest point on one of the lines
        # (added to the file as fields by the Near operation above
        new_x = row.GetValue("NEAR_X")
        new_y = row.GetValue("NEAR_Y")

        # Create a new point object with the new x and y values
        point = gp.CreateObject("Point")
        point.x = new_x
        point.y = new_y

        # Assign it to the shape field
        row.shape = point

        # Update the row data and move to the next row
        rows.UpdateRow(row)
        row = rows.Next()

def snap_points_alt(points, lines, distance):
    """
    Ogi's updated snap_points function.
    """

    points = arcpy.Near_analysis(points, lines, str(distance), "LOCATION")

    # Create an update cursor for the points Feature Class
    # making sure that the NEAR_X and NEAR_Y fields are included
    # in the return data
    with arcpy.da.UpdateCursor(points, ["NEAR_X", "NEAR_Y", "SHAPE@XY"]) as cursor:
        for row in cursor:
            x, y, shape_xy = row
            shape_xy = (x, y)
            cursor.updateRow([x, y, shape_xy])
    return(points)

### Local variables:
## input

# Necessary to change
env.workspace = 'C:/ecan/local/Projects/Waimakariri/GIS/vector'

boundary = 'C:/ecan/local/Projects/Waimakariri/GIS/vector/waimak_cwms.shp'

# May not be necessary to change
recorders = 'S:/Surface Water/shared\\GIS_base\\vector\\rec_nat_sites1.shp'
min_sites = 'S:/Surface Water/shared/GIS_base/vector/low_flows/gauge_nat_sites1.shp'
streams = 'S:/Surface Water/shared\\GIS_base\\vector\\MFE_REC_rivers_no_1st.shp'
dem = 'S:/Surface Water/shared\\GIS_base\\raster\\DEM_8m_2012\\linz_8m_dem'

point_dis = 1000
stream_depth = 10
grid_size = 8
pour_dis = 20
env.extent = boundary

merge_cond = "ReferenceN \"ReferenceN\" true true false 12 Text 0 0 ,First,#,min_gauge_sites_otop,ReferenceN,-1,-1;SITENUMBER \"SITENUMBER\" true true false 10 Long 0 10 ,First,#,otop_recorders_v02,SITENUMBER,-1,-1"

## output

min_sites_loc = 'gauge_nat_sites1.shp'
recorders_loc = 'rec_nat_sites1.shp'
streams_loc = 'MFE_streams_loc.shp'
dem_loc = 'dem_loc.tif'
min_gauge_sites = 'gauge_nat_sites1.shp'
sites = 'sites2.shp'
stream_diss = 'MFE_rivers_diss.shp'
stream_rast = 'stream_rast2.tif'
accu_tif = 'accu1.tif'
catch_poly = 'catch2.shp'
catch_sites_join = 'catch_sites_join.shp'
catch_sites_csv = 'results/catch_sites.csv'
sites_csv = 'results/sites.csv'
catch_poly_csv = 'results/catch.csv'

##########################
#### Processing

### Process stie and streams vectors

# Clip lots of things
arcpy.Clip_analysis(recorders, boundary, recorders_loc)
arcpy.Clip_analysis(min_sites, boundary, min_sites_loc)
arcpy.Clip_analysis(streams, boundary, streams_loc)


# Select the gauging sites from the min flow sites
arcpy.Select_analysis(min_sites_loc, min_gauge_sites, "ReferenceS = 'Gauging'")

# Merge min_gauge_sites to recorders
arcpy.Merge_management([min_gauge_sites, recorders_loc], sites, merge_cond)

# create site numeric field
arcpy.AddField_management(sites, "site_num", "LONG")
arcpy.CalculateField_management(sites, "site_num", '!SITENUMBER! if !SITENUMBER! > 0 else int(!ReferenceN!)', "PYTHON_9.3")

# Snap sites to streams
snap_points(sites, streams_loc, point_dis)

# Process: Dissolve
arcpy.Dissolve_management(streams_loc, stream_diss, "", "", "MULTI_PART", "DISSOLVE_LINES")

# Add raster parameters to streams layer
arcpy.AddField_management(stream_diss, "rast", "SHORT")
arcpy.CalculateField_management(stream_diss, "rast", stream_depth, "PYTHON_9.3")

############################################
### Delineate catchments

# Clip the DEM to the study area
arcpy.Clip_management(dem, "1323813.1799 5004764.9257 1688157.0305 5360238.95", dem_loc, boundary, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# Convert stream vector to raster
arcpy.FeatureToRaster_conversion(stream_diss, 'rast', stream_rast, grid_size)

# Fill holes in DEM
dem_fill = Fill(dem_loc)

# Subtract stream raster from
s_rast = Raster(stream_rast)
dem_diff = Con(IsNull(s_rast), dem_fill, dem_fill - s_rast)

# Fill holes in DEM
dem2 = Fill(dem_diff)

# flow direction
fd1 = FlowDirection(dem2)

# flow accu
accu1 = FlowAccumulation(fd1)
#ah.FlowAccumulation(fd1, accu2)
accu1.save(accu_tif)
accu1 = 'accu1.tif'

# create pour points
pp1 = SnapPourPoint(sites, accu1, pour_dis, 'site')

# Determine the catchments for all points
catch1 = Watershed(fd1, pp1)

# Convert raster to polygon
arcpy.RasterToPolygon_conversion(catch1, catch_poly, 'SIMPLIFY', 'Value')

# Add in a field for the area of each catchment
arcpy.AddField_management(catch_poly, "area_m2", "LONG")
arcpy.CalculateField_management(catch_poly, "area_m2", 'round(!shape.area!)', "PYTHON_9.3")

## Remove tiny polygons
#arcpy.SelectLayerByAttribute_management(catch_poly, "NEW_SELECTION", '"area_m2 < 2000')

### STOP and check that the sites have a realistic catchment area!

##############################################
### Spatial join to determine which site is upstream of each catchment area

arcpy.SpatialJoin_analysis(catch_poly, sites, catch_sites_join, "JOIN_ONE_TO_MANY", "KEEP_ALL", "GRIDCODE \"GRIDCODE\" true true false 10 Long 0 10 ,First,#,catch1,GRIDCODE,-1,-1;site \"site\" true true false 10 Long 0 10 ,First,#,site,site,-1,-1", "WITHIN_A_DISTANCE", str(pour_dis + 10) + " Meters", "")

############################################
#### Export data
arcpy.ExportXYv_stats(catch_sites_join, "GRIDCODE;site", "COMMA", catch_sites_csv, "ADD_FIELD_NAMES")
arcpy.ExportXYv_stats(sites, "SITENUMBER;ReferenceN;site", "COMMA", sites_csv, "ADD_FIELD_NAMES")
arcpy.ExportXYv_stats(catch_poly, "ID;GRIDCODE;area_m2", "COMMA", catch_poly_csv, "ADD_FIELD_NAMES")


