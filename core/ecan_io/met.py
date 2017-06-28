# -*- coding: utf-8 -*-
"""
Functions for importing meteorological data.
"""


def rd_vcn(select='all', data_dir='Y:/VirtualClimate/VCN_precip_ET_2016-06-06', data_type='precip', site_col='ID', comp_table='Y:/VirtualClimate/VCN_id_comp_table.csv', id_col=1, buffer_dis=0, vcsn_grid_shp='Y:/VirtualClimate/NIWA_rain_grid_Canterbury.shp', vcsn_site_col='Data VCN_s', netcdf_export=False, csv_export=False, export_path='vcn_data.csv'):
    """
    Function to read many VCN csv files and concatenate them into a dataframe where each column is a VCN station.

    Arguments:\n
    data_dir -- The directory with all of the VCN csv files.\n
    select -- Either 'all', a path to a csv that contains a 'site_col' with the station numbers that should be included, or a list of the site numbers.\n
    site_col -- If 'select' is not 'all', then the column name of the station numbers for the sites csv.\n
    data_type -- Either 'precip' or 'ET', not both.\n
    comp_table -- The path to the table that has both station ID types.\n
    export -- Should the data be exported?\n
    export_path -- The path where the data should be exported to.
    """
    from os import path
    from numpy import in1d, ndarray, tile, repeat
    from core.misc import rd_dir
    from pandas import read_csv, concat, to_numeric, melt
    from geopandas import read_file
    import xarray as xr

    spatial_bool = False
    files, sites = rd_dir(data_dir, 'csv', True)

    if (type(select) is str) and (select is not 'all'):
        if select.endswith('.shp'):
            spatial_bool = True

            #### Read in data
            poly1 = read_file(select)
            if type(id_col) is str:
                poly2 = poly1[[id_col, 'geometry']]
            else:
                id_col_name = poly1.columns[id_col - 1]
                poly2 = poly1[[id_col_name, 'geometry']]
            poly2.columns = ['id', 'geometry']

            vcn_grid = read_file(vcsn_grid_shp)[[vcsn_site_col, 'geometry']]
            vcn_grid.columns = ['station', 'geometry']

            #### Perform vector operations for initial processing
            ## Dissolve polygons by id
            poly3 = poly2.dissolve(by='id')
            poly0 = poly3.unary_union

            ## Create buffer
            poly_buff = poly0.buffer(buffer_dis)

            ## Select only the vcn sites within the buffer
            vcn_grid2 = vcn_grid[vcn_grid.within(poly_buff)]
            site_index = in1d(sites, vcn_grid2.station.values)
            vcn_grid3 = vcn_grid2[in1d(vcn_grid2.station.values, sites)].to_crs(epsg=4326)
            files = files[site_index]
            sites1 = sites[site_index]
            t1 = read_csv(comp_table)
            sites = t1.loc[in1d(t1.ecan_id, sites1)].sort_values('ecan_id').net_id.values

            ## Extract the x/y coordinates
            x = vcn_grid3.geometry.apply(lambda p: p.x).round(3).values
            y = vcn_grid3.geometry.apply(lambda p: p.y).round(3).values

        else:
            sites1 = read_csv(select)[site_col]
            first1 = sites1[0]
            if 'P' in first1:
                t1 = read_csv(comp_table)
                t2 = t1.loc[in1d(t1.net_id, sites1), 'ecan_id']
                index_sites = in1d(sites, t2.values)
            else:
                index_sites = in1d(sites, sites1.astype('int').values)
            sites = sites1.values
            files = files[index_sites]
    elif (type(select) is list) or (type(select) is ndarray):
        first1 = select[0]
        if 'P' in first1:
            t1 = read_csv(comp_table)
            t2 = t1.loc[in1d(t1.net_id, select), 'ecan_id']
            index_sites = in1d(sites, t2.values)
        else:
            index_sites = in1d(sites, select)
        sites = select
        files = files[index_sites]

    if netcdf_export:
        df1 = concat((read_csv(path.join(data_dir, f), usecols=['precip', 'ET']) for f in files), axis=0)
        time_index1 = read_csv(path.join(data_dir, files[0]), index_col=0, parse_dates=True, infer_datetime_format=True).index
        time_index = tile(time_index1, len(files))
        site_index = repeat(sites, len(time_index1))

        df1.set_index([time_index, site_index], inplace=True)
        df1.index.names = ['time', 'site']
        df2 = concat([to_numeric(df1[d], errors='coerce') for d in df1.columns], axis=1)

        et_df = df2['ET'].reset_index().pivot('time', 'site')
        precip_df = df2['precip'].reset_index().pivot('time', 'site')

        xr1 = xr.Dataset({'precip': (['time', 'site'], precip_df.values), 'ET': (['time', 'site'], et_df.values)}, coords={'site': sites, 'time': time_index1.values, 'x': ('site', x), 'y': ('site', y)})

        xr1.to_netcdf(export_path)

        if spatial_bool:
            vcn_grid3 = vcn_grid2.geometry
            vcn_grid3.index = vcn_grid2.station
            return([xr1, vcn_grid3])
        else:
            return(xr1)
    else:
        df1 = concat((read_csv(path.join(data_dir, f), index_col=0, parse_dates=True, infer_datetime_format=True)[data_type] for f in files), axis=1)
        df1.columns = sites
        if csv_export:
            df1.to_csv(export_path)
        if spatial_bool:
            vcn_grid3 = vcn_grid2.geometry
            vcn_grid3.index = vcn_grid2.station
            return([df1, vcn_grid3])
        else:
            return(df1)


def proc_metservice_nc(nc, lat_coord='south_north', lon_coord='west_east', time_coord='Time', time_var='Times'):
    """
    Function to process MetService netcdf files so that it is actually complete. The function adds in the appropriate coordinate arrays for the data and resaves the file with '_corr" added to the end of the name.

    nc -- Full path to the MetService nc file (str).\n
    lat_coord -- The name of the lat coordinate that should be added (str).\n
    lon_coord -- Same as lat_coord except for the lon.\n
    time_coord -- Ditto for the time.\n
    time_var -- The existing name of the time variable (that should be converted and removed).
    """
    from xarray import open_dataset
    from os import path
    from numpy import arange
    from pandas import to_datetime
    from core.ecan_io.met import ACPR_to_rate
    from core.spatial.vector import convert_crs

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
#
#    ### Put in the hourly rate
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
    new_path = path.splitext(nc)[0] + '_corr.nc'
    x5.to_netcdf(new_path)
    x1.close()
    x5.close()


def ACPR_to_rate(df, lat_coord='y', lon_coord='x', time_coord='time'):
    """
    Function to convert cummulative precip to hourly rate.

    df -- DataFrame of the cummulative precip.\n
    lat_coord -- The name of the lat coordinate that should be added (str).\n
    lon_coord -- Same as lat_coord except for the lon.\n
    time_coord -- Ditto for the time.
    """
    from pandas import merge

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

    nc -- The path to the corrected MetService netcdf file.\n
    lat_coord -- The name of the lat coordinate that should be added (str).\n
    lon_coord -- Same as lat_coord except for the lon.\n
    time_coord -- Ditto for the time.\n
    precip_var -- The precip variable name.\n
    proj4 -- The proj4 coordinate system attribute name.
    """
    from xarray import open_dataset
    from numpy import tile
    from shapely.geometry import Point
    from geopandas import GeoDataFrame

    ### Extract all data to dataframes
    ds = open_dataset(nc)
    precip = ds[precip_var].to_dataframe().reset_index()
    proj1 = str(ds.attrs[proj4])

    ### Create geodataframe
    time = precip[time_coord].unique()
    sites0 = precip.loc[precip[time_coord] == time[0], [lon_coord, lat_coord]]
    precip.loc[:, 'site'] = tile(sites0.index, len(time))
    sites0.index.name = 'site'

    geometry = [Point(xy) for xy in zip(sites0[lon_coord], sites0[lat_coord])]
    sites = GeoDataFrame(sites0.index, geometry=geometry, crs=proj1)

    ### Return
    ds.close()
    return(precip, sites)


def sel_interp_agg(precip, sites, poly, grid_res, data_col, time_col, x_col, y_col, buffer_dis=10000, interp_fun='multiquadric', agg_ts_fun=None, period=None, digits=3, agg_xy=False, output_format='csv', nfiles='many', output_path='precip_interp.csv'):
    """
    Function to select the precip sites within a polygon with a certain buffer distance, then interpolate/resample the data at a specific resolution, then output the results.
    precip -- dataframe of time, x, y, and precip.\n
    sites -- GeoDataFrame of site locations.\n
    poly -- String path of a shapefile polygon.\n
    res -- Resolution in meters of the resampling.\n
    buffer_dis -- Buffer distance of the polygon selection.\n
    interp_fun -- The scipy Rbf interpolation function to be applied (see https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.Rbf.html).\n
    agg_ts_fun -- The pandas time series resampling function to resample the data in time (either 'mean' or 'sum'). If None, then no time resampling.\n
    agg_ts_fun -- The pandas time series code to resample the data in time (i.e. '2H' for two hours).\n
    digits -- the number of digits to round to (int).\n
    agg_xy -- Should all of the interpolated points within the polygon area be aggregated (mean) to a single time series?\n
    output_format -- Either a str or list of 'csv', 'geotiff', and/or 'netcdf'.\n
    nfiles -- If 'geotiff' is in the output_format, then 'one' or 'many' geotiffs to be created.\n
    output_path -- Full path string where the output should be stored.
    """

    from core.spatial import sel_sites_poly, grid_interp_ts, xy_to_gpd, save_geotiff
    from geopandas import read_file
    from numpy import tile
    from os import path

    ### Select the locations within the polygon
    poly1 = read_file(poly)
    sites1 = sites.to_crs(poly1.crs)
    sites_sel = sel_sites_poly(sites1, poly, buffer_dis)
    sites2 = sites[sites.site.isin(sites_sel.site)]

    ### Select the precip data from the sites
    precip2 = precip[precip.site.isin(sites2.site)]

    ### Interpolate grid
    poly_crs = ['+' + str(i) + '=' + str(poly1.crs[i]) for i in poly1.crs]
    poly_crs1 = ' '.join(poly_crs)
    new_precip = grid_interp_ts(precip2, time_col, x_col, y_col, data_col, grid_res, sites.crs, poly_crs1, interp_fun=interp_fun, agg_ts_fun=agg_ts_fun, period=period, digits=digits)

    ### Create new sites list
    time = new_precip[time_col].sort_values().unique()
    sites_new_df = new_precip.loc[new_precip[time_col] == time[0], [x_col, y_col, data_col]]
    sites_new = xy_to_gpd(sites_new_df.index.values, x_col, y_col, sites_new_df, poly_crs1)
    sites_new.columns = ['site', 'geometry']
    new_precip['site'] = tile(sites_new_df.index.values, len(time))

    ### Select sites from polygon
    sites_sel2 = sel_sites_poly(sites_new, poly)
    new_precip2 = new_precip.loc[new_precip.site.isin(sites_sel2.site), [time_col, x_col, y_col, data_col]]

    ### Agg to polygon if required
    if agg_xy:
        new_precip3 = new_precip2.groupby(time_col)[data_col].mean().round(digits)
        time_col = None
    else:
        new_precip3 = new_precip2.set_index([time_col, x_col, y_col])[data_col]

    ### Save results
    path1 = path.splitext(output_path)[0]
    if 'csv' in output_format:
        new_precip3.to_csv(path1 + '.csv', header=True)

    if 'geotiff' in output_format:
        df = new_precip3.reset_index()
        save_geotiff(df=df, data_col=data_col, crs=poly_crs1, x_col=x_col, y_col=y_col, time_col=time_col, nfiles=nfiles, export_path=path1 + '.tif')

    if 'netcdf' in output_format:
        ds1 = new_precip3.to_xarray().to_dataset()
        ds1.attrs['spatial_ref'] = poly_crs1
        ds1.to_netcdf(path1 + '.nc')

    return(new_precip3)


def proc_niwa_rcp(base_path, mtypes, poly, vcsn_sites_csv=r'Z:\Data\VirtualClimate\GIS\niwa_vcsn_wgs84.csv', id_col='Network', x_col='deg_x', y_col='deg_y', output_fun=None, export_path='output'):
    """
    Function to read in the NIWA RCP netcdf files and output the data in a specified format.
    """
    from pandas import read_csv
    from core.spatial import xy_to_gpd, sel_sites_poly
    from geopandas import read_file
    from os import path, walk, makedirs
    from core.ecan_io.met import rd_niwa_rcp_dir

    mtype_name = {'precip': 'TotalPrecipCorr', 'T_max': 'MaxTempCorr', 'T_min': 'MinTempCorr', 'P_atmos': 'MSLP', 'PET': 'PE', 'RH_mean': 'RelHum', 'R_s': 'SurfRad', 'U_z': 'WindSpeed'}

    ### Import and reorganize data
    vcsn_sites = read_csv(vcsn_sites_csv)[[id_col, x_col, y_col]]

    sites_gpd = xy_to_gpd(id_col, x_col, y_col, vcsn_sites, 4326)
    poly1 = read_file(poly)

    sites_gpd2 = sites_gpd.to_crs(poly1.crs)

    mtypes1 = [mtype_name[i] for i in mtypes]

    ### Select sites
    sites_gpd3 = sel_sites_poly(sites_gpd2, poly1)[id_col]
    site_loc1 = vcsn_sites[vcsn_sites[id_col].isin(sites_gpd3)]
    site_loc1.columns = ['id', 'x', 'y']

    ### Read and extract data from netcdf files

    for root, dirs, files in walk(base_path):
        files2 = [i for i in files if i.endswith('.nc')]
        files3 = [j for j in files2 if any(j.startswith(i) for i in mtypes1)]
        file_paths1 = [path.join(root, i) for i in files3]
        if len(file_paths1) > 0:
            ds = rd_niwa_rcp_dir(file_paths1, site_loc1, mtypes)
            if callable(output_fun):
                new_base_path = root.replace(base_path, export_path)
                base_file_name = file_paths1[0].split('VCSN_')[1]
                if not path.exists(new_base_path):
                    makedirs(new_base_path)
                output_fun(ds, new_base_path, base_file_name)
                print(base_file_name)
            else:
                raise ValueError('Must have a output function.')

    ### What should I return?


def rd_niwa_rcp_dir(file_paths, site_loc, mtypes):
    """
    Function to read in one or more nc files with the same time, x, and y but different mtypes.

    file_paths -- A string of a file path or a list of string paths.\n
    site_loc -- A dataframe with id and x and y in decimal degrees WGS84.\n
    mtypes -- The measurement types to extract.
    """
    from xarray import open_dataset, Dataset, DataArray
    from numpy import in1d
    from os.path import basename

    ### Parameters
    mtype_param = {'precip': 'rain', 'T_max': 'tmax', 'T_min': 'tmin', 'P_atmos': 'mslp', 'PET': 'pe', 'RH_mean': 'rh', 'R_s': 'srad', 'U_z': 'wind'}
    mtype_param1 = {v: k for k, v in mtype_param.iteritems()}
    prob_mtypes = ['P_atmos', 'RH_mean', 'R_s', 'U_z', 'T_min']
    bad_names = {'mslp2': 'mslp'}
    data_attr = {'grid_mapping': 'crs'}
    nc_crs = {'inverse_flattening': 298.257223563, 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137, 'grid_mapping_name': 'latitude_longitude'}
    time_attr = {'bounds': 'time_bounds', 'standard_name': 'time', 'axis': 'T', 'long_name': 'time (end of reporting period)'}

    ### Extract the proper time coordinate to fix the problem parameters if needed
    if any(in1d(prob_mtypes, mtypes)):
        tmax_file = [j for j in file_paths if 'MaxTempCorr' in basename(j)][0]
        tmax_ds = open_dataset(tmax_file)
        time_da = tmax_ds['time'].copy()
        tmax_ds.close()

    ### Open data files and Extract the data
    ds10 = Dataset()
    for i in file_paths:
        ds5 = open_dataset(i)
        if any(in1d(bad_names.keys(), ds5.data_vars.keys())):
            ds5 = ds5.rename(bad_names)
        mtype2 = [j for j in ds5.data_vars.keys() if j in mtype_param.values()][0]
        mtype0 = mtype_param1[mtype2]

        ## Prepare the selection from the x and y
        if 'bool_lat' not in locals():
            lat1 = (ds5.latitude.data * 1000).astype('int32')
            lon1 = (ds5.longitude.data * 1000).astype('int32')
            site_lat  = (site_loc['y'] * 1000).astype('int32').unique()
            site_lon  = (site_loc['x'] * 1000).astype('int32').unique()

            bool_lat = in1d(lat1, site_lat)
            bool_lon = in1d(lon1, site_lon)

        ## Extract the data based on criteria from earlier
        ds6 = ds5.sel(latitude=bool_lat, longitude=bool_lon)
        da1 = ds6[[mtype2]]
        attr1 = da1[mtype2].attrs
        attr1.update(data_attr)
        da1[mtype2].attrs = attr1
        da1['time'].attrs = time_attr

        ## Imbed the correct time coordinates if necessary
        if mtype0 in prob_mtypes:
            da1['time'] = time_da

        ## Merge datasets
        ds10 = ds10.merge(da1)
#        print([mtype2, len(ds10.time)])
        print(mtype0)

    ### Add in dummy GIS variable
    ds_crs = DataArray(4326, attrs=nc_crs, name='crs').to_dataset()
    ds11 = ds10.merge(ds_crs).copy()
    ds11.attrs = da1.attrs

    return(ds11)


#def export_rcp_lst(ds, export_path):
#    """
#    Function to take the output of rd_niwa_rcp_dir and save the data as standard lst files.
#    """
#    from os import path
#
#    ### Reorganize
#    df3 = df[['id', 'y', 'x', 'time', 'precip', 'PET']]
#    time1 = df3.time.dt.strftime('%Y%m%d')
#    df3.loc[:, 'time'] = time1
#
#    ### Save to many files (by id)
#    id1 = df3.id.unique()
#    for i in id1:
#        out1 = df3[df3.id == i]
#        out1.to_csv(path.join(export_path, i + '.lst'), header=False, index=False)


def export_rcp_nc(ds, export_path, file_name):
    """
    Function to take the output of rd_niwa_rcp_dir and save the data as a standard nc file.
    """
    from os.path import join

    ### Save to nc file based on directory names
    ds.to_netcdf(join(export_path, 'VCSN_' + file_name))
    ds.close()



def nc_add_gis(nc, x_coord, y_coord):
    """
    Function to add the appropriate attributes to a netcdf file to be able to load it into GIS if the netcdf file has x and y in WGS84 decimal degrees.

    nc -- A path str to the netcdf file (str).\n
    x_coord -- The x coordinate name (str).\n
    y_coord -- The y coordinate name (str).
    """
    from xarray import open_dataset, DataArray
    from os.path import splitext

    ### Attributes for the various datasets
    nc_crs = {'inverse_flattening': 298.257223563, 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137, 'grid_mapping_name': 'latitude_longitude'}

    x_attr = {'long_name': 'longitude', 'units': 'degrees_east', 'standard_name': 'longitude', 'axis': 'X'}
    y_attr = {'long_name': 'latitude', 'units': 'degrees_north', 'standard_name': 'latitude', 'axis': 'Y'}
    data_attr = {'grid_mapping': 'crs'}

    ### Read in the nc
    ds1 = open_dataset(nc)

    ### Determine the variables with x and y coordinates
    vars1 = ds1.data_vars
    vars2 = [i for i in vars1 if ((x_coord in ds1[i]) & (y_coord in ds1[i]))]

    ### Put in the additional attribute into the variables
    ds1[x_coord].attrs = x_attr
    ds1[y_coord].attrs = y_attr

    for i in vars2:
        attr1 = ds1[i].attrs
        attr1.update(data_attr)
        ds1[i].attrs = attr1

    ### Add crs dummy dataset
    ds_crs = DataArray(4326, attrs=nc_crs, name='crs').to_dataset()
    ds2 = ds1.merge(ds_crs)

    ### Resave nc file
    new_path = splitext(nc)[0] + '_gis.nc'
    ds2.to_netcdf(new_path)
    ds1.close()
    ds2.close()


def sel_niwa_vcsn(mtypes, sites, nc_path='Y:/VirtualClimate/vcsn_precip_et_2016-06-06.nc', vcsn_sites_csv=r'Y:\VirtualClimate\GIS\niwa_vcsn_wgs84.csv', id_col='Network', x_col='deg_x', y_col='deg_y', include_sites=False, out_crs=None, netcdf_out=None):
    """
    Function to read in the NIWA vcsn netcdf file and output the data as a dataframe.

    mtypes -- A string or list of the measurement types (either 'precip', or 'PET').\n
    sites -- Either a list of vcsn site names or a polygon of the area of interest.\n
    nc_path -- The path to the vcsn nc file.\n
    vcsn_sites_csv -- The csv file that relates the site name to coordinates.\n
    id_col -- The site name column in vcsn_sites_csv.\n
    x_col - The x column name in vcsn_sites_csv.\n
    y_col -- The y column name in vcsn_sites_csv.\n
    include_sites -- Should the site names be added to the output?\n
    out_crs -- The crs epsg number for the output coordinates if different than the default WGS85 (e.g. 2193 for NZTM).
    """
    from pandas import read_csv, Series, merge
    from core.spatial import xy_to_gpd, sel_sites_poly, convert_crs
    from geopandas import read_file
    from numpy import ndarray
    from xarray import open_dataset

    mtype_name = {'precip': 'rain', 'PET': 'pe'}

    ### Import and reorganize data
    vcsn_sites = read_csv(vcsn_sites_csv)[[id_col, x_col, y_col]]

    if isinstance(sites, str):
        if sites.endswith('.shp'):
            sites_gpd = xy_to_gpd(id_col, x_col, y_col, vcsn_sites, 4326)
            poly1 = read_file(sites)

            sites_gpd2 = sites_gpd.to_crs(poly1.crs)

            ### Select sites
            sites2 = sel_sites_poly(sites_gpd2, poly1)[id_col]
    elif isinstance(sites, (list, Series, ndarray)):
        sites2 = sites

    ### Select locations
    site_loc1 = vcsn_sites[vcsn_sites[id_col].isin(sites2)]
    site_loc1.columns = ['id', 'x', 'y']

    ### Select mtypes
    if isinstance(mtypes, str):
        mtypes1 = [mtype_name[mtypes]]
    else:
        mtypes1 = [mtype_name[i] for i in mtypes]

    if include_sites:
        mtypes1.extend(['site'])

    ### Read and extract data from netcdf files
    ds1 = open_dataset(nc_path)
    ds2 = ds1.sel(longitude=site_loc1.x.unique(), latitude=site_loc1.y.unique())
    ds3 = ds2[mtypes1]

    ### Convert to DataFrame
    df1 = ds3.to_dataframe().reset_index()
    df1.rename(columns={'latitude': 'y', 'longitude': 'x'}, inplace=True)

    ### Convert to different crs if needed
    if out_crs is not None:
        crs1 = convert_crs(out_crs)
        new_gpd1 = xy_to_gpd('id', 'x', 'y', site_loc1, 4326)
        new_gpd2 = new_gpd1.to_crs(crs1)
        site_loc2 = site_loc1.copy()
        site_loc2['x_new'] = new_gpd2.geometry.apply(lambda j: j.x)
        site_loc2['y_new'] = new_gpd2.geometry.apply(lambda j: j.y)

        df2 = merge(df1, site_loc2[['x', 'y', 'x_new', 'y_new']], on=['x', 'y'])
        df3 = df2.drop(['x', 'y'], axis=1).rename(columns={'x_new': 'x', 'y_new': 'y'})
        col_order = ['y', 'x', 'time']
        col_order.extend(mtypes1)
        df4 = df3[col_order]
    else:
        df4 = df1

    ds1.close()
    ds3.close()

    ### Return
    if isinstance(netcdf_out, str):
        ds3.to_netcdf(netcdf_out)
    return(df4)















