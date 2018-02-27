# -*- coding: utf-8 -*-
"""
Functions for importing data from various sources.
"""
import xarray as xr
from core.spatial import xy_to_gpd, pts_poly_join
from geopandas import read_file
from pandas import merge


#####################################
### General data import


def rd_nc(poly_shp, nc_path, poly_epsg=4326, poly_id='Station_ID', x_col='longitude', y_col='latitude', data_col='rain',
          as_ts=True, export=True, export_path='nc_data.csv'):
    """
    Function to read in netCDF files, select locations based on a polygon, and export the results.
    """

    ### Read in all data
    poly = read_file(poly_shp)[[poly_id, 'geometry']].to_crs(epsg=poly_epsg)
    nc = xr.open_dataset(nc_path)

    ### Filter nc data
    df1 = nc.to_dataframe().drop('time_bnds', axis=1).reset_index()
    df1 = df1[df1.nb2 == 0].drop('nb2', axis=1)

    ### convert x and y to geopandas
    df1_xy = df1[[y_col, x_col]].drop_duplicates()
    df1_xy['id'] = range(len(df1_xy))
    pts = xy_to_gpd('id', x_col, y_col, df1_xy, poly_epsg)

    ### Mask the points from the polygon
    join1, poly2 = pts_poly_join(pts, poly, poly_id)
    join2 = join1[['id', poly_id]]

    ### Select the associated data
    sel_xy = merge(df1_xy, join2, on='id').drop('id', axis=1)
    df2 = merge(df1, sel_xy, on=[y_col, x_col])

    ### Convert to time series
    if as_ts:
        df3 = df2[[poly_id, 'time', data_col]].groupby([poly_id, 'time']).first().reset_index()
        df4 = df3.pivot(index='time', columns=poly_id, values=data_col).round(2)
        if export:
            df4.to_csv(export_path)
    else:
        df4 = df2
        if export:
            df4.to_csv(export_path)

    return df4
