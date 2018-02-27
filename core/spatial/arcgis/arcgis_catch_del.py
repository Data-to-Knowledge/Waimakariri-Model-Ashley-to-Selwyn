# -*- coding: utf-8 -*-
"""
Functions to delineate catchments using arcpy.
"""


def arc_catch_del(WD, boundary_shp, sites_shp, site_num_col='site', point_dis=1000, stream_depth=10, grid_size=8, pour_dis=20, streams='S:/Surface Water/shared\\GIS_base\\vector\\MFE_REC_rivers_no_1st.shp', dem='S:/Surface Water/shared\\GIS_base\\raster\\DEM_8m_2012\\linz_8m_dem', export_dir='results', overwrite_rasters=False):
    """
    Arcpy function to delineate catchments based on specific points, a polygon, and the REC rivers layer.

    Be careful that the folder path isn't too long!!! Do not have spaces in the path name!!! Arc sucks!!!

    WD -- Working directory (str).\n
    boundary_shp -- The path to the shapefile polygon boundary extent (str).\n
    sites_shp -- The path to the sites shapefile (str).\n
    site_num_col -- The column in the sites_shp that contains the site IDs (str).\n
    point_dis -- The max distance to snap the sites to the nearest stream line (int).\n
    stream_depth -- The depth that the streams shapefile should be burned into the dem (int).\n
    grid_size -- The resolution of the dem (int).\n
    streams -- The path to the streams shapefile (str).\n
    dem -- The path to the dem (str).\n
    export_dir -- The subfolder where the results should be saved (str).\n
    overwrite_rasters -- Should the flow direction and flow accumulation rasters be overwritten? (bool).
    """

    # load in the necessary arcpy libraries to import arcpy
    import sys
    sys.path.append('C:\\Python27\\ArcGIS10.4\\Lib\\site-packages')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\arcpy')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\ArcToolbox\Scripts')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\bin')
    sys.path.append('C:\\Python27\\ArcGIS10.4\\lib')

    # Import packages
    import arcpy
    from arcpy import env
    from arcpy.sa import Raster, Con, IsNull, FlowDirection, FlowAccumulation, Fill, SnapPourPoint, Watershed
    import os
    #import ArcHydroTools as ah

    # Check out spatial analyst license
    arcpy.CheckOutExtension('Spatial')
    # Define functions

#    def snap_points(points, lines, distance):
#
#        import arcgisscripting, sys
#
#        gp = arcgisscripting.create()
#
#        # Load the Analysis toolbox so that the Near tool is available
#        gp.toolbox = "analysis"
#
#        # Perform the Near operation looking for the nearest line
#        # (from the lines Feature Class) to each point (from the
#        # points Feature Class). The third argument is the search
#        # radius - blank means to search as far as is needed. The
#        # fourth argument instructs the command to output the
#        # X and Y co-ordinates of the nearest point found to the
#        # NEAR_X and NEAR_Y fields of the points Feature Class
#        gp.near(points, lines, str(distance), "LOCATION")
#
#        # Create an update cursor for the points Feature Class
#        # making sure that the NEAR_X and NEAR_Y fields are included
#        # in the return data
#        rows = gp.UpdateCursor(points, "", "", "NEAR_X, NEAR_Y")
#
#        row = rows.Next()
#
#        # For each row
#        while row:
#            # Get the location of the nearest point on one of the lines
#            # (added to the file as fields by the Near operation above
#            new_x = row.GetValue("NEAR_X")
#            new_y = row.GetValue("NEAR_Y")
#
#            # Create a new point object with the new x and y values
#            point = gp.CreateObject("Point")
#            point.x = new_x
#            point.y = new_y
#
#            # Assign it to the shape field
#            row.shape = point
#
#            # Update the row data and move to the next row
#            rows.UpdateRow(row)
#            row = rows.Next()

    def snap_points(points, lines, distance):
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

    ### Parameters:
    ## input

    # Necessary to change
    env.workspace = WD
    boundary = boundary_shp
    sites_in = sites_shp

#    site_num_col = 'site'

    # May not be necessary to change
    final_export_dir = export_dir
#    streams = 'S:/Surface Water/shared\\GIS_base\\vector\\MFE_REC_rivers_no_1st.shp'
#    dem = 'S:/Surface Water/shared\\GIS_base\\raster\\DEM_8m_2012\\linz_8m_dem'

    env.extent = boundary
    arcpy.env.overwriteOutput = True

    ## output
    bound = 'bound_diss.shp'
    sites = 'sites_bound.shp'
    streams_loc = 'MFE_streams_loc.shp'
    dem_loc = 'dem_loc.tif'
    stream_diss = 'MFE_rivers_diss.shp'
    stream_rast = 'stream_rast.tif'
    dem_diff_tif = 'dem_diff.tif'
    dem_fill_tif = 'dem_fill.tif'
    fd_tif = 'fd1.tif'
    accu_tif = 'accu1.tif'
    catch_poly = 'catch_del.shp'

    if not os.path.exists(os.path.join(env.workspace, final_export_dir)):
        os.makedirs(os.path.join(env.workspace, final_export_dir))

    ##########################
    #### Processing

    ### Process sites and streams vectors

    # Dissolve boundary for faster processing
    arcpy.Dissolve_management(boundary, bound)

    # Clip sites and streams to boundary
    arcpy.Clip_analysis(streams, bound, streams_loc)
    arcpy.Clip_analysis(sites_in, bound, sites)

    # Snap sites to streams layer
    snap_points(sites, streams_loc, point_dis)

    # Dissolve stream network
    arcpy.Dissolve_management(streams_loc, stream_diss, "", "", "MULTI_PART", "DISSOLVE_LINES")

    # Add raster parameters to streams layer
    arcpy.AddField_management(stream_diss, "rast", "SHORT")
    arcpy.CalculateField_management(stream_diss, "rast", stream_depth, "PYTHON_9.3")

    ############################################
    ### Delineate catchments

    # Convert stream vector to raster
    arcpy.FeatureToRaster_conversion(stream_diss, 'rast', stream_rast, grid_size)

    ## Create the necessary flow direction and accumulation rasters if they do not already exist
    if os.path.exists(os.path.join(env.workspace, accu_tif)) & (not overwrite_rasters):
        accu1 = Raster(accu_tif)
        fd1 = Raster(fd_tif)
    else:
        # Clip the DEM to the study area
        print('clipping DEM to catchment area...')
        arcpy.Clip_management(dem, "1323813.1799 5004764.9257 1688157.0305 5360238.95", dem_loc, bound, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

        # Fill holes in DEM
        print('Filling DEM...')
#        dem_fill = Fill(dem_loc)

        # Subtract stream raster from
        s_rast = Raster(stream_rast)
        dem_diff = Con(IsNull(s_rast), dem_loc, dem_loc - s_rast)
        dem_diff.save(dem_diff_tif)

        # Fill holes in DEM
        dem2 = Fill(dem_diff_tif)
        dem2.save(dem_fill_tif)

        # flow direction
        print('Flow direction...')
        fd1 = FlowDirection(dem2)
        fd1.save(fd_tif)

        # flow accu
        print('Flow accumulation...')
        accu1 = FlowAccumulation(fd1)
        accu1.save(accu_tif)

    # create pour points
    pp1 = SnapPourPoint(sites, accu1, pour_dis, site_num_col)

    # Determine the catchments for all points
    catch1 = Watershed(fd1, pp1)

    # Convert raster to polygon
    arcpy.RasterToPolygon_conversion(catch1, catch_poly, 'SIMPLIFY', 'Value')

    # Add in a field for the area of each catchment
    arcpy.AddField_management(catch_poly, "area_m2", "LONG")
    arcpy.CalculateField_management(catch_poly, "area_m2", 'round(!shape.area!)', "PYTHON_9.3")

    #### Check back in the spatial analyst license once done
    arcpy.CheckInExtension('Spatial')


##############################################
### Spatial join to determine which site is upstream of each catchment area


def arc_spatial_join(WD, site_num_col='site', pour_dis=20, export_dir='results', catch_sites_csv = 'catch_sites.csv', catch_poly_csv = 'catch.csv'):
    # load in the necessary arcpy libraries to import arcpy
    import sys
    sys.path.append('C:\\Python27\\ArcGIS10.4\\Lib\\site-packages')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\arcpy')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\ArcToolbox\Scripts')
    sys.path.append(r'C:\Program Files (x86)\ArcGIS\Desktop10.4\bin')
    sys.path.append('C:\\Python27\\ArcGIS10.4\\lib')

    # Import packages
    import arcpy
    from arcpy import env
    import os

    env.workspace = WD
    final_export_dir = export_dir

    catch_poly = 'catch_del.shp'
    sites = 'sites_bound.shp'
    catch_sites_join = 'catch_sites_join.shp'

    arcpy.SpatialJoin_analysis(catch_poly, sites, catch_sites_join, "JOIN_ONE_TO_MANY", "KEEP_ALL", "", "WITHIN_A_DISTANCE", str(pour_dis + 10) + " Meters", "")

    # Remove unnecessary fields
    keep_fields = ['FID', 'Shape', 'GRIDCODE', site_num_col]
    rem_fields = [f.name for f in arcpy.ListFields(catch_sites_join)]
    [rem_fields.remove(x) for x in keep_fields]

    arcpy.DeleteField_management(catch_sites_join, rem_fields)

    ############################################
    #### Export data
    arcpy.ExportXYv_stats(catch_sites_join, "GRIDCODE;" + site_num_col, "COMMA", os.path.join(final_export_dir, catch_sites_csv), "ADD_FIELD_NAMES")
    arcpy.ExportXYv_stats(catch_poly, "ID;GRIDCODE;area_m2", "COMMA", os.path.join(final_export_dir, catch_poly_csv), "ADD_FIELD_NAMES")

    ###########################################
    #### Check back in the spatial analyst license once done
    arcpy.CheckInExtension('Spatial')



