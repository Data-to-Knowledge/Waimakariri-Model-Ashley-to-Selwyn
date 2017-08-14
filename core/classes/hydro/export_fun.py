# -*- coding: utf-8 -*-
"""
Functions for exporting data from within a hydro class.
"""
from numpy import in1d
from pandas import DataFrame, concat

def to_csv(self, csv_path, mtypes=None, sites=None, pivot=False, resample=True, require=None):
    """
    Function to export data from a hydro class to a MultiIndex csv.
    """

    df1 = self.sel_ts(mtypes=mtypes, sites=sites, pivot=pivot, resample=resample, require=require)
    df1.to_csv(csv_path, header=True)


def to_netcdf(self, nc_path):
    """
    Function to export a copy of a hydro class object to a netcdf file.
    """
    from xarray import Dataset, DataArray

    ### package ts data
    ds0 = Dataset(self.data.reset_index(['site', 'time']))

    ### package site attribue data
    if hasattr(self, 'site_attr'):
        site_attr = getattr(self, 'site_attr').copy()
        site_attr.columns = ['site_attr_' + i for i in site_attr.columns]
        site_attr.index.name = 'site_attr_site'
        ds0 = ds0.merge(Dataset(site_attr))

    ### package geo data
    if hasattr(self, 'geo_loc'):
        geo1 = getattr(self, 'geo_loc').copy()
        geo1['x'] = geo1.geometry.apply(lambda x: x.x)
        geo1['y'] = geo1.geometry.apply(lambda x: x.y)
        geo2 = geo1.drop('geometry', axis=1)
        geo2.columns = ['geo_loc_' + i for i in geo2.columns]
        geo2.index.name = 'geo_loc_' + geo2.index.name
        ds0 = ds0.merge(Dataset(geo2))
        ds0.attrs['crs'] = geo2.crs

    if hasattr(self, 'geo_catch'):
        geo1 = getattr(self, 'geo_catch').copy()
        geo1['wkt'] = geo1.geometry.apply(lambda x: x.to_wkt())
        geo2 = geo1.drop('geometry', axis=1)
        geo2.columns = ['geo_catch_' + i for i in geo2.columns]
        geo2.index.name = 'geo_catch_' + geo2.index.name
        ds0 = ds0.merge(Dataset(geo2))

    ### Export data
    ds0.to_netcdf(nc_path)
    ds0.close()


def to_shp(self, shp_path):
    """
    Function to export a shapefile of the site locations.
    """

    if not hasattr(self, 'geo_loc'):
        raise ValueError('Object has no geo locations!')

    ### Prepare output
    geo1 = self.geo_loc
    sites = self.sites
    mtypes_sites = self.mtypes_sites

    if len(sites) != len(geo1):
        geo1.to_file(shp_path)
    else:
        t1 = {i: in1d(sites, list(mtypes_sites[i])).astype(str).tolist() for i in mtypes_sites}
        df1 = DataFrame(t1, index=sites)
        geo2 = concat([geo1, df1], axis=1)
        geo2.reset_index().to_file(shp_path)













