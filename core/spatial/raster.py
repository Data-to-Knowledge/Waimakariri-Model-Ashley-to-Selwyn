# -*- coding: utf-8 -*-
"""
Raster and spatial interpolation functions.
"""
from pandas import DataFrame, TimeGrouper, Grouper, to_datetime
from numpy import tile, repeat, column_stack, nan, arange, meshgrid
from shapely.geometry import Point
from geopandas import GeoDataFrame, read_file
from scipy.interpolate import griddata, Rbf
from core.spatial.vector import convert_crs
from rasterio import open as ras_open
from rasterio import transform
from os import path


def grid_interp_ts(df, time_col, x_col, y_col, data_col, grid_res, from_crs=None, to_crs=2193, interp_fun='cubic', agg_ts_fun=None, period=None, digits=2):
    """
    Function to take a dataframe of z values and interate through and resample both in time and space. Returns a DataFrame structured like df.

    df -- DataFrame containing four columns as shown in the below parameters.\n
    time_col -- The time column name.\n
    x_col -- The x column name.\n
    y_col -- The y column name.\n
    data_col -- The data column name.\n
    grid_res -- The resulting grid resolution in meters (or the unit of the final projection).\n
    from_crs -- The projection info for the input data if the result should be reprojected to the to_crs projection (either a proj4 str or epsg int).\n
    to_crs -- The projection for the output data similar to from_crs.\n
    interp_fun -- The scipy Rbf interpolation function to be applied (see https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.Rbf.html).\n
    agg_ts_fun -- The pandas time series resampling function to resample the data in time (either 'mean' or 'sum'). If None, then no time resampling.\n
    period -- The pandas time series code to resample the data in time (i.e. '2H' for two hours).\n
    digits -- the number of digits to round to (int).
    """

    #### Create the grids
    df1 = df.copy()

    #### Resample the time series data
    if agg_ts_fun is not None:
        df1a = df1.set_index(time_col)
        if agg_ts_fun == 'sum':
            df2 = df1a.groupby([TimeGrouper(period), Grouper(y_col), Grouper(x_col)])[data_col].sum().reset_index()
        elif agg_ts_fun == 'mean':
            df2 = df1a.groupby([TimeGrouper(period), Grouper(y_col), Grouper(x_col)])[data_col].mean().reset_index()
        else:
            raise ValueError("agg_ts_fun should be either 'sum' or 'mean'.")
        time = df2[time_col].unique()
    else:
        df2 = df1

    time = df2[time_col].sort_values().unique()

    if from_crs is None:
        x = df2.loc[df2[time_col] == time[0], x_col].values
        y = df2.loc[df2[time_col] == time[0], y_col].values
    else:
        data1 = df2.loc[df2[time_col] == time[0]]
        from_crs1 = convert_crs(from_crs, pass_str=True)
        to_crs1 = convert_crs(to_crs, pass_str=True)
        geometry = [Point(xy) for xy in zip(data1[x_col], data1[y_col])]
        gpd = GeoDataFrame(data1.index, geometry=geometry, crs=from_crs1)
        gpd1 = gpd.to_crs(crs=to_crs1)
        x = gpd1.geometry.apply(lambda p: p.x).round(digits).values
        y = gpd1.geometry.apply(lambda p: p.y).round(digits).values

    xy = column_stack((x, y))

    max_x = x.max()
    min_x = x.min()

    max_y = y.max()
    min_y = y.min()

    new_x = arange(min_x, max_x, grid_res)
    new_y = arange(min_y, max_y, grid_res)
    x_int, y_int = meshgrid(new_x, new_y)

    #### Create new df
    x_int2 = x_int.flatten()
    y_int2 = y_int.flatten()
    xy_int = column_stack((x_int2, y_int2))
    time_df = repeat(time, len(x_int2))
    x_df = tile(x_int2, len(time))
    y_df = tile(y_int2, len(time))
    new_df = DataFrame({'time': time_df, 'x': x_df, 'y': y_df, data_col: repeat(0, len(time) * len(x_int2))})

    new_lst = []
    for t in to_datetime(time):
        set1 = df2.loc[df2[time_col] == t, data_col]
#        index = new_df[new_df['time'] == t].index
        new_z = griddata(xy, set1.values, xy_int, method=interp_fun).round(digits)
        new_z[new_z < 0] = 0
        new_lst.extend(new_z.tolist())
#        print(t)
    new_df.loc[:, data_col] = new_lst

    #### Export results
    return(new_df[new_df[data_col].notnull()])


def point_interp_ts(df, time_col, x_col, y_col, data_col, point_shp, point_site_col, from_crs, to_crs=None, interp_fun='cubic', agg_ts_fun=None, period=None, digits=2):
    """
    Function to take a dataframe of z values and interate through and resample both in time and space. Returns a DataFrame structured like df.

    df -- DataFrame containing four columns as shown in the below parameters.\n
    time_col -- The time column name (str).\n
    x_col -- The x column name (str).\n
    y_col -- The y column name (str).\n
    data_col -- The data column name (str).\n
    from_crs -- The projection info for the input data (either a proj4 str or epsg int).\n
    point_shp -- Path to shapefile of points to be interpolated (str) or a GeoPandas GeoDataFrame.\n
    point_site_col -- The column name of the site names/numbers of the point_shp (str).\n
    to_crs -- The projection for the output data similar to from_crs.\n
    interp_fun -- The scipy Rbf interpolation function to be applied (see https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.Rbf.html).\n
    agg_ts_fun -- The pandas time series resampling function to resample the data in time (either 'mean' or 'sum'). If None, then no time resampling.\n
    period -- The pandas time series code to resample the data in time (i.e. '2H' for two hours).\n
    digits -- the number of digits to round to (int).
    """

    #### Read in points
    if isinstance(point_shp, str) & isinstance(point_site_col, str):
        points = read_file(point_shp)[[point_site_col, 'geometry']]
        to_crs1 = points.crs
    elif isinstance(point_shp, GeoDataFrame) & isinstance(point_site_col, str):
        points = point_shp[[point_site_col, 'geometry']]
        to_crs1 = points.crs
    else:
        raise ValueError('point_shp must be a str path to a shapefile or a GeoDataFrame and point_site_col must be a str.')

    #### Create the grids
    df1 = df.copy()

    #### Resample the time series data
    if agg_ts_fun is not None:
        df1a = df1.set_index(time_col)
        if agg_ts_fun == 'sum':
            df2 = df1a.groupby([TimeGrouper(period), Grouper(y_col), Grouper(x_col)])[data_col].sum().reset_index()
        elif agg_ts_fun == 'mean':
            df2 = df1a.groupby([TimeGrouper(period), Grouper(y_col), Grouper(x_col)])[data_col].mean().reset_index()
        else:
            raise ValueError("agg_ts_fun should be either 'sum' or 'mean'.")
        time = df2[time_col].unique()
    else:
        df2 = df1

    time = df2[time_col].sort_values().unique()

    #### Convert input data to crs of points shp and create input xy
    data1 = df2.loc[df2[time_col] == time[0]]
    from_crs1 = convert_crs(from_crs, pass_str=True)

    if to_crs is not None:
        to_crs1 = convert_crs(to_crs, pass_str=True)
        points = points.to_crs(to_crs1)
    geometry = [Point(xy) for xy in zip(data1[x_col], data1[y_col])]
    gpd = GeoDataFrame(data1.index, geometry=geometry, crs=from_crs1)
    gpd1 = gpd.to_crs(crs=to_crs1)
    x = gpd1.geometry.apply(lambda p: p.x).round(digits).values
    y = gpd1.geometry.apply(lambda p: p.y).round(digits).values

    xy = column_stack((x, y))

    #### Prepare the x and y of the points geodataframe output
    x_int = points.geometry.apply(lambda p: p.x).round(digits).values
    y_int = points.geometry.apply(lambda p: p.y).round(digits).values
    sites = points[point_site_col]

    xy_int = column_stack((x_int, y_int))

    #### Create new df
    sites_ar = tile(sites, len(time))
    time_ar = repeat(time, len(xy_int))
    x_ar = tile(x_int, len(time))
    y_ar = tile(y_int, len(time))
    new_df = DataFrame({'site': sites_ar, 'time': time_ar, 'x': x_ar, 'y': y_ar, data_col: repeat(0, len(time) * len(xy_int))})

    new_lst = []
    for t in to_datetime(time):
        set1 = df2.loc[df2[time_col] == t, data_col]
        new_z = griddata(xy, set1.values, xy_int, method=interp_fun).round(digits)
        new_z[new_z < 0] = 0
        new_lst.extend(new_z.tolist())
#        print(t)
    new_df.loc[:, data_col] = new_lst

    #### Export results
    return(new_df[new_df[data_col].notnull()])


def grid_resample(x, y, z, x_int, y_int, digits=3, method='multiquadric'):
    """
    Function to interpolate and resample a set of x, y, z values.
    """

    interp1 = Rbf(x, y, z, function=method)
    z_int = interp1(x_int, y_int).round(digits)
    z_int[z_int < 0] = 0

    z_int2 = z_int.flatten()
    return(z_int2)


def save_geotiff(df, data_col, crs, x_col='x', y_col='y', time_col=None, nfiles='many', export_path='geotiff.tif', grid_res=None):
    """
    Function to convert a dataframe of x, y, and data to a GeoTiff. If the DataFrame has a time_col, then these instances can be saved as multiple bands in the GeoTiff or as multiple GeoTiffs.

    df -- DataFrame with at least an x_col, y_col, and data_col. The combo of x_col and y_col must be unique without a corresponding time_col.\n
    x_col -- The x column.\n
    y_col -- The y column.\n
    data_col - The data column.\n
    crs -- The projection info for the data (either a proj4 str or epsg int).\n
    time_col -- The time column if one exists.\n
    nfiles -- If time_col is passed, how many files should be created? 'one' will make a single GeoTiff with many bands and 'many' will make many GeoTiffs.\n
    path -- The save path.\n
    grid_res -- The grid resolution of the output raster (int). The default None will output the the resolution based on the point spacing of a regular grid.
    """

    ### create the xy coordinates
    if time_col is None:
        xy1 = df[[x_col, y_col]]
    else:
        time = df[time_col].sort_values().unique()
        xy1 = df.loc[df[time_col] == time[0], [x_col, y_col]]
    if any(xy1.duplicated()):
        raise ValueError('x and y coordinates are not unique!')

    ### Determine grid res
    if grid_res is None:
        res_df1 = (xy1.loc[0] - xy1).abs()
        res_df2 = res_df1.replace(0, nan).min()
        x_res = res_df2[x_col]
        y_res = res_df2[y_col]
    elif isinstance(grid_res, int):
        x_res = y_res = grid_res
    else:
        raise ValueError('grid_res must either be None or an integer.')

    ### Make the affline transformation for Rasterio
    trans2 = transform.from_origin(xy1[x_col].min() - x_res/2, xy1[y_col].max() + y_res/2, x_res, y_res)

    ### Make the rasters
    if time_col is None:
        z = df.set_index([y_col, x_col])[data_col].unstack().values[::-1]
        new_dataset = ras_open(export_path, 'w', driver='GTiff', height=len(xy1[y_col].unique()), width=len(xy1[x_col].unique()), count=1, dtype=df[data_col].dtype, crs=convert_crs(crs, pass_str=True), transform=trans2)
        new_dataset.write(z, 1)
        new_dataset.close()
    else:
        if nfiles == 'one':
            new_dataset = ras_open(export_path, 'w', driver='GTiff', height=len(xy1[y_col].unique()), width=len(xy1[x_col].unique()), count=len(time), dtype=df[data_col].dtype, crs=convert_crs(crs), transform=trans2)
            for i in range(1, len(time)+1):
                z = df.loc[df[time_col] == time[i - 1]].set_index([y_col, x_col])[data_col].unstack().values[::-1]
                new_dataset.write(z, i)
            new_dataset.close()
        elif nfiles == 'many':
            file1 = path.splitext(export_path)[0]
            for i in time:
                str_date = to_datetime(i).strftime('%Y-%m-%d_%H')
                file2 = file1 + '_' + str_date + '.tif'
                z = df.loc[df[time_col] == i].set_index([y_col, x_col])[data_col].unstack().values[::-1]
                new_dataset = ras_open(file2, 'w', driver='GTiff', height=len(xy1[y_col].unique()), width=len(xy1[x_col].unique()), count=1, dtype=df[data_col].dtype, crs=convert_crs(crs), transform=trans2)
                new_dataset.write(z, 1)
                new_dataset.close()






