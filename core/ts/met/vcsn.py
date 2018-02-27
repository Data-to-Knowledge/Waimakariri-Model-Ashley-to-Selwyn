# -*- coding: utf-8 -*-
"""
Functions for processing vcsn met data from NIWA.
"""
from xarray import open_dataset, DataArray
from os.path import splitext


def nc_add_gis(nc, x_coord, y_coord):
    """
    Function to add the appropriate attributes to a netcdf file to be able to load it into GIS if the netcdf file has x and y in WGS84 decimal degrees.

    nc -- A path str to the netcdf file (str).\n
    x_coord -- The x coordinate name (str).\n
    y_coord -- The y coordinate name (str).
    """

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




