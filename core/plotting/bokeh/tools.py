# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 09:39:40 2017

@author: MichaelEK
"""
from shapely.geometry import Point, LineString, Polygon


def getPolyCoords(row, coord_type, geom='geometry'):
    """Returns the coordinates ('x' or 'y') of edges of a Polygon exterior"""

    # Parse the exterior of the coordinate
    exterior = row[geom].exterior

    if coord_type == 'x':
        # Get the x coordinates of the exterior
        return list(exterior.coords.xy[0])
    elif coord_type == 'y':
        # Get the y coordinates of the exterior
        return list(exterior.coords.xy[1])


def getPointCoords(row, coord_type, geom='geometry'):
    """Calculates coordinates ('x' or 'y') of a Point geometry"""
    if coord_type == 'x':
        return row[geom].x
    elif coord_type == 'y':
        return row[geom].y


def getLineCoords(row, coord_type, geom='geometry'):
    """Returns a list of coordinates ('x' or 'y') of a LineString geometry"""
    if coord_type == 'x':
        return list(row[geom].coords.xy[0])
    elif coord_type == 'y':
        return list(row[geom].coords.xy[1])


def getCoords(gdf):
    """
    Function to put in x and y coordinates for Bokeh plotting from geodataframes of points, lines, or polygons.
    """
    gdf1 = gdf.copy()

    x_all = []
    y_all = []
    for i in gdf1.index:
        geo1 = gdf1.geometry[i]
        if isinstance(geo1, Point):
            x = geo1.x
            y = geo1.y
        elif isinstance(geo1, LineString):
            x = list(geo1.coords.xy[0])
            y = list(geo1.coords.xy[1])
        elif isinstance(geo1, Polygon):
            x = list(geo1.exterior.coords.xy[0])
            y = list(geo1.exterior.coords.xy[1])
        else:
            raise TypeError('Index ' + str(i) + ' is not a shapely Point, LineString, or Polygon')
        x_all.append(x)
        y_all.append(y)

    gdf1.loc[:, 'x'] = x_all
    gdf1.loc[:, 'y'] = y_all

    return(gdf1)





















