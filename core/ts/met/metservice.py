# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:17:26 2017

@author: MichaelEK

Functions for processing MetService data.
"""
from os import path
from xarray import open_dataset
from numpy import tile, arange
from shapely.geometry import Point
from geopandas import GeoDataFrame
from pandas import to_datetime, merge, to_numeric
from core.ecan_io.mssql import rd_sql
from core.spatial.vector import xy_to_gpd


def proc_metservice_nc(nc, lat_coord='south_north', lon_coord='west_east', time_coord='Time', time_var='Times', export_dir=None):
    """
    Function to process MetService netcdf files so that it is actually complete. The function adds in the appropriate coordinate arrays for the data and resaves the file with '_corr" added to the end of the name.

    nc : str
        Full path to the MetService nc file (str).
    lat_coord : str
        The name of the lat coordinate that should be added (str).
    lon_coord : str
        Same as lat_coord except for the lon.
    time_coord : str
        Ditto for the time.
    time_var : str
        The existing name of the time variable (that should be converted and removed).
    export_dir : str or None
        The export directory for the processed netcdf file. If None, then the new file is put into the same directory as the original file.

    Returns
    -------
    str
        The new path to the processed netcdf file.
    """

    ### Parameters
    proj1 = '+proj=lcc +lat_1=-60 +lat_2=-30 +lat_0=-60 +lon_0=167.5 +x_0=211921 +y_0=-1221320 +a=6367470 +b=6367470 +no_defs'

    ### Read in the nc file
    x1 = open_dataset(nc)

    ### Extract parameters and convert to numpy arrays
    time1 = to_datetime(x1[time_var].data, format='%Y-%m-%d_%H:%M:%S')

    nlat = x1.dims[lat_coord]
    nlon = x1.dims[lon_coord]
    x_res = int(x1.attrs['DX'])
    y_res = int(x1.attrs['DY'])

    lat = arange(nlat, dtype='int32') * y_res
    lon = arange(nlon, dtype='int32') * x_res

    ### Remove the old time variable and add in the coordinates
    x2 = x1.drop(time_var)
    x2.coords[time_coord] = ((time_coord), time1)
    x2.coords[lat_coord] = ((lat_coord), lat)
    x2.coords[lon_coord] = ((lon_coord), lon)

    ### rename coordinates
    x3 = x2.rename({time_coord: 'time', lat_coord: 'y', lon_coord: 'x'})

    ### Calc hourly precip rate
    df = x3['ACPR'].to_dataframe().reset_index()
    precip = ACPR_to_rate(df, 'y', 'x')

    ### Remove the first time step (as there is no data for it)
    x4 = x3.sel(time=precip.time.unique())

    ### Put in the hourly rate
    precip_ds = precip.set_index(['time', 'y', 'x']).to_xarray()
    x5 = x4.merge(precip_ds)

    ### Add in attributes
    ## x
    x_attrs = {'standard_name': 'projection_x_coordinate', 'units': 'm', 'axis': 'X'}
    x5.coords['x'].attrs = x_attrs

    ## y
    y_attrs = {'standard_name': 'projection_y_coordinate', 'units': 'm', 'axis': 'Y'}
    x5.coords['y'].attrs = y_attrs

    ## variables
    ACPR_attrs = {'standard_name': 'precipitation_amount', 'units': 'mm', 'description': 'accumulated total grid precipitation'}
    precip_attrs = {'standard_name': 'precipitation_amount', 'units': 'mm', 'description': 'hourly precipitation'}

    x5.variables['ACPR'].attrs = ACPR_attrs
    x5.variables['precip'].attrs = precip_attrs

    ## Overall attributes
    x5.attrs['spatial_ref'] =  proj1

    ### Save the new file and close them
    if export_dir is None:
        new_path = path.splitext(nc)[0] + '_corr.nc'
    elif isinstance(export_dir, (str, unicode)):
        nc_file = path.splitext(path.split(nc)[1])[0] + '_corr.nc'
        new_path = path.join(export_dir, nc_file)

    x5.to_netcdf(new_path)
    x1.close()
    x5.close()
    return(new_path)


def ACPR_to_rate(df, lat_coord='y', lon_coord='x', time_coord='time'):
    """
    Function to convert cummulative precip to hourly rate.

    df : DataFrame
        DataFrame of the cummulative precip.
    lat_coord : str
        The name of the lat coordinate that should be added (str).
    lon_coord : str
        Same as lat_coord except for the lon.
    time_coord : str
        Ditto for the time.

    Returns
    -------
    DataFrame
        Three dimensions with hourly precip rate.
    """

    ### Extract data into dataframe
    df1 = df.copy().set_index(time_coord)
    df1a = df1.shift(1, freq='H')
    df0 = merge(df1.reset_index(), df1a.reset_index(), on=[time_coord, lon_coord, lat_coord], how='inner')
    df0['precip'] = (df0['ACPR_x'] - df0['ACPR_y']).round(3)
    df2 = df0[[time_coord, lon_coord, lat_coord, 'precip']]

    return(df2)


def MetS_nc_to_df(nc, lat_coord='y', lon_coord='x', time_coord='time', precip_var='precip', proj4='spatial_ref'):
    """
    Function to convert a MetService nc file to the components of precip and sites with x y locations.

    nc : str
        The path to the corrected MetService netcdf file.
    lat_coord : str
        The name of the lat coordinate that should be added (str).
    lon_coord : str
        Same as lat_coord except for the lon.
    time_coord : str
        Ditto for the time.
    precip_var : str
        The precip variable name.
    proj4 : str
        The proj4 coordinate system attribute name.

    Returns
    -------
    DataFrame
        Precip rate by time, x, and y (with site)
    GeoDataFrame
        Site locations dataframe
    Timestamp
        The model date start time
    """

    ### Extract all data to dataframes
    with open_dataset(nc) as ds:
        precip = ds[precip_var].to_dataframe().reset_index()
    proj1 = str(ds.attrs[proj4])

    ### Create geodataframe
    time = precip[time_coord].unique()
    sites0 = precip.loc[precip[time_coord] == time[0], [lon_coord, lat_coord]]
    precip.loc[:, 'site'] = tile(sites0.index, len(time))
    sites0.index.name = 'site'

    geometry = [Point(xy) for xy in zip(sites0[lon_coord], sites0[lat_coord])]
    sites = GeoDataFrame(sites0.index, geometry=geometry, crs=proj1)

    start_date = to_datetime(ds.attrs['START_DATE'], format='%Y-%m-%d_%H:%M:%S')

    ### Return
    ds.close()
    return(precip, sites, start_date)


def metconnect_id_loc(sites=None, mc_server='SQL2012PROD03', mc_db='MetConnect', mc_site_table='RainFallPredictionSites', mc_cols=['MetConnectID', 'SiteString', 'TidedaID'], gis_server='SQL2012PROD05'):
    """
    Function to extract the metconnect id table with geometry location.

    Parameters
    ----------
    sites : list of int or None
        The site numbers to extract from the table, or None for all.

    Returns
    -------
    GeoDataFrame
    """

    ### Input parameters
#    hy_server = 'SQL2012PROD05'
#    hy_db = 'Hydrotel'
#    pts_table = 'Points'
#    objs_table = 'Objects'
#    sites_table = 'Sites'
#
#    pts_cols = ['Point', 'Object']
#    objs_cols = ['Object', 'Site']
#    sites_cols = ['Site', 'ExtSysId']

    loc_db = 'Bgauging'
    loc_table = 'RSITES'

    loc_cols = ['SiteNumber', 'NZTMX', 'NZTMY']

    ## Import tables
    mc1 = rd_sql(mc_server, mc_db, mc_site_table, mc_cols)
    mc2 = mc1[~mc1.SiteString.str.startswith('M')]
    mc2.columns = ['MetConnectID', 'site_name', 'ExtSysId']
    mc2 = mc2[(mc2.MetConnectID != 7) & mc2.ExtSysId.notnull()]
    mc2.loc[:, 'ExtSysId'] = mc2.loc[:, 'ExtSysId'].astype(int)

#    hy_pts = rd_sql(hy_server, hy_db, pts_table, pts_cols, 'Point', mc2.Point.tolist())
#    hy_objs = rd_sql(hy_server, hy_db, objs_table, objs_cols, 'Object', hy_pts.Object.tolist())
#    hy_sites = rd_sql(hy_server, hy_db, sites_table, sites_cols, 'Site', hy_objs.Site.tolist())
#    hy_sites['ExtSysId'] = to_numeric(hy_sites['ExtSysId'])
    hy_loc = rd_sql(gis_server, loc_db, loc_table, loc_cols, 'SiteNumber', mc2.ExtSysId.tolist())
    hy_loc.columns = ['ExtSysId', 'x', 'y']

#    t1 = merge(mc2, hy_pts, on='Point')
#    t2 = merge(t1, hy_objs, on='Object')
#    t3 = merge(t2, hy_sites, on='Site')
    t4 = merge(mc2, hy_loc, on='ExtSysId')

    hy_xy = xy_to_gpd('MetConnectID', 'x', 'y', t4)

    return(hy_xy)























