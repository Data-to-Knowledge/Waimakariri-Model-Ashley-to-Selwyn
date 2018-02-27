# -*- coding: utf-8 -*-
"""
Functions for importing meteorological data.
"""


def rd_niwa_rcp(base_path, mtypes, poly,
                vcsn_sites_csv=r'\\fileservices02\ManagedShares\Data\VirtualClimate\GIS\niwa_vcsn_wgs84.csv',
                id_col='Network', x_col='deg_x', y_col='deg_y', output_fun=None, export_path='output'):
    """
    Function to read in the NIWA RCP netcdf files and output the data in a specified format.
    """
    from pandas import read_csv
    from core.spatial import xy_to_gpd, sel_sites_poly
    from geopandas import read_file
    from os import path, walk, makedirs
    from core.ecan_io.met import rd_niwa_rcp_dir

    mtype_name = {'precip': 'TotalPrecipCorr', 'T_max': 'MaxTempCorr', 'T_min': 'MinTempCorr', 'P_atmos': 'MSLP',
                  'PET': 'PE', 'RH_mean': 'RelHum', 'R_s': 'SurfRad', 'U_z': 'WindSpeed'}

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
    mtype_param = {'precip': 'rain', 'T_max': 'tmax', 'T_min': 'tmin', 'P_atmos': 'mslp', 'PET': 'pe', 'RH_mean': 'rh',
                   'R_s': 'srad', 'U_z': 'wind'}
    mtype_param1 = {v: k for k, v in mtype_param.iteritems()}
    prob_mtypes = ['P_atmos', 'RH_mean', 'R_s', 'U_z', 'T_min']
    bad_names = {'mslp2': 'mslp'}
    data_attr = {'grid_mapping': 'crs'}
    nc_crs = {'inverse_flattening': 298.257223563, 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137,
              'grid_mapping_name': 'latitude_longitude'}
    time_attr = {'bounds': 'time_bounds', 'standard_name': 'time', 'axis': 'T',
                 'long_name': 'time (end of reporting period)'}

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
            site_lat = (site_loc['y'] * 1000).astype('int32').unique()
            site_lon = (site_loc['x'] * 1000).astype('int32').unique()

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

    return (ds11)


# def export_rcp_lst(ds, export_path):
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
    nc_crs = {'inverse_flattening': 298.257223563, 'longitude_of_prime_meridian': 0, 'semi_major_axis': 6378137,
              'grid_mapping_name': 'latitude_longitude'}

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


def rd_niwa_vcsn(mtypes, sites,
                 nc_path=r'\\fileservices02\ManagedShares\Data\VirtualClimate\vcsn_precip_et_2016-06-06.nc',
                 vcsn_sites_csv=r'\\fileservices02\ManagedShares\Data\VirtualClimate\GIS\niwa_vcsn_wgs84.csv',
                 id_col='Network', x_col='deg_x', y_col='deg_y', buffer_dis=0, include_sites=False, from_date=None,
                 to_date=None, out_crs=None, netcdf_out=None):
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
    from pandas import read_csv, Series, merge, to_datetime
    from core.spatial import xy_to_gpd, sel_sites_poly, convert_crs
    from geopandas import read_file
    from numpy import ndarray, in1d
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
            sites2 = sel_sites_poly(sites_gpd2, poly1, buffer_dis)[id_col]
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
    time1 = to_datetime(ds1.time.values)
    if isinstance(from_date, str):
        time1 = time1[time1 >= from_date]
    if isinstance(to_date, str):
        time1 = time1[time1 <= to_date]
    lat1 = ds1.latitude.values
    lon1 = ds1.longitude.values
    lat2 = lat1[in1d(lat1, site_loc1.y.unique())]
    lon2 = lon1[in1d(lon1, site_loc1.x.unique())]
    ds2 = ds1.loc[{'longitude': lon2, 'time': time1.values, 'latitude': lat2}]
    ds3 = ds2[mtypes1]

    ### Convert to DataFrame
    df1 = ds3.to_dataframe().reset_index()
    df1.rename(columns={'latitude': 'y', 'longitude': 'x'}, inplace=True)
    df1 = df1.dropna()

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
    return (df4)


def rd_niwa_climate_proj(mtypes, bound_shp, nc_dir, buffer_dis=0, from_date=None, to_date=None, out_crs=None,
                         netcdf_out=None):
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
    from core.misc import rd_dir
    from os.path import join
    from core.ecan_io.met import sel_xy_nc

    file_name = {'precip': 'TotalPrecipCorr', 'PET': 'PE'}
    mtype_name = {'precip': 'rain', 'PET': 'pe'}

    ### Select mtypes
    if isinstance(mtypes, str):
        mtypes1 = [mtypes]
    else:
        mtypes1 = mtypes
    niwa_mtypes = [mtype_name[i] for i in mtypes1]

    ### Read and extract data from netcdf files
    files1 = rd_dir(nc_dir, 'nc')
    for mt in mtypes1:
        niwa_mtype = mtype_name[mt]
        niwa_file_name = file_name[mt]
        file1 = [i for i in files1 if niwa_file_name in i][0]

        ds3 = sel_xy_nc(bound_shp, join(nc_dir, file1), nc_vars=[niwa_mtype], buffer_dis=buffer_dis,
                        from_date=from_date, to_date=to_date, out_type='xarray')

        if 'ds4' in locals():
            ds4 = ds4.merge(ds3)
        else:
            ds4 = ds3

    ### Convert to DataFrame
    df1 = ds4.to_dataframe().reset_index()
    df1.rename(columns={'latitude': 'y', 'longitude': 'x'}, inplace=True)
    df1 = df1.dropna()
    df1.set_index('time', inplace=True)
    df1_grp = df1.groupby(['x', 'y'])
    df2 = df1_grp.resample('D')[niwa_mtypes].first().interpolate().round(1).reset_index()

    ### Return
    if isinstance(netcdf_out, str):
        ds4.to_netcdf(netcdf_out)
    ds4.close()
    return (df2)


def rd_niwa_data_lsrm(bound_shp, nc_dir,
                      vcsn_sites_csv=r'\\fileservices02\ManagedShares\Data\VirtualClimate\GIS\niwa_vcsn_wgs84.csv',
                      id_col='Network', x_col='deg_x', y_col='deg_y', buffer_dis=0, from_date=None, to_date=None,
                      out_crs=None, netcdf_out=None):
    """
    Convienence function to read in either the past VCSN data or the projections.
    """
    from core.ecan_io.met import rd_niwa_climate_proj, rd_niwa_vcsn
    from os import path

    ### Determine whether the input is past data or projections
    past_vcsn = 'vcsn_precip_et_2016-06-06.nc'
    name_split = path.split(nc_dir)[1]

    ### Run function to get data
    if name_split == past_vcsn:
        ds5 = rd_niwa_vcsn(['precip', 'PET'], bound_shp, buffer_dis=buffer_dis, from_date=from_date, to_date=to_date,
                           out_crs=out_crs, netcdf_out=netcdf_out)
    else:
        ds5 = rd_niwa_climate_proj(['precip', 'PET'], bound_shp, nc_dir=nc_dir, buffer_dis=buffer_dis,
                                   from_date=from_date, to_date=to_date, out_crs=out_crs, netcdf_out=netcdf_out)
    return (ds5)


def sel_xy_nc(bound_shp, nc_path, x_col='longitude', y_col='latitude', time_col='time', nc_vars=None, buffer_dis=0,
              from_date=None, to_date=None, nc_crs=4326, out_crs=None, out_type='pandas'):
    """
    Function to select space and time data from a netcdf file using a polygon shapefile.
    """
    from pandas import Series, merge, to_datetime
    from core.spatial import xy_to_gpd, convert_crs
    from geopandas import read_file
    from numpy import ndarray
    from xarray import open_dataset

    ### Process the boundary layer
    bound = read_file(bound_shp).buffer(buffer_dis).to_crs(convert_crs(nc_crs))
    x_min, y_min, x_max, y_max = bound.unary_union.bounds

    ### Read and extract data from netcdf files
    ds1 = open_dataset(nc_path)
    time1 = to_datetime(ds1[time_col].values)
    if isinstance(from_date, str):
        time1 = time1[time1 >= from_date]
    if isinstance(to_date, str):
        time1 = time1[time1 <= to_date]
    lat1 = ds1[y_col].values
    lon1 = ds1[x_col].values
    lat2 = lat1[(lat1 >= y_min) & (lat1 <= y_max)]
    lon2 = lon1[(lon1 >= x_min) & (lon1 <= x_max)]
    ds2 = ds1.loc[{x_col: lon2, time_col: time1.values, y_col: lat2}]

    #    coords1 = ds2.coords.keys()
    #    dims1 = ds2.dims.keys()

    ## Select mtypes
    if isinstance(nc_vars, str):
        ds3 = ds2[[nc_vars]]
    elif isinstance(nc_vars, (list, ndarray, Series)):
        ds3 = ds2[nc_vars]
    elif nc_vars is None:
        ds3 = ds2

    ### Convert to different crs if needed
    if out_crs is not None:
        df1 = ds3.to_dataframe().reset_index()
        xy1 = ds3[[x_col, y_col]].copy()
        xy2 = xy1.to_dataframe().reset_index()
        crs1 = convert_crs(out_crs)
        new_gpd1 = xy_to_gpd(xy2.index, x_col, y_col, xy2, nc_crs)
        new_gpd2 = new_gpd1.to_crs(crs1)
        site_loc2 = xy2.copy()
        site_loc2['x_new'] = new_gpd2.geometry.apply(lambda j: j.x)
        site_loc2['y_new'] = new_gpd2.geometry.apply(lambda j: j.y)

        df2 = merge(df1, site_loc2[[x_col, y_col, 'x_new', 'y_new']], on=[x_col, y_col], how='left')
        df3 = df2.drop([x_col, y_col], axis=1).rename(columns={'x_new': x_col, 'y_new': y_col})
        ds1.close()
        return (df3)
    elif out_type == 'pandas':
        df1 = ds3.to_dataframe().reset_index()
        ds1.close()
        return (df1)
    elif out_type == 'xarray':
        return (ds3)
