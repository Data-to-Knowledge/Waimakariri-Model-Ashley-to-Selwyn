# -*- coding: utf-8 -*-
"""
Functions for exporting data from within a hydro class.
"""
from numpy import in1d
from pandas import DataFrame, concat
from xarray import Dataset, DataArray
from core.spatial.vector import convert_crs

def to_csv(self, csv_path, mtypes=None, sites=None, pivot=False, resample=None, require=None):
    """
    Function to export data from a hydro class to a MultiIndex csv.
    """

    df1 = self.sel_ts(mtypes=mtypes, sites=sites, pivot=pivot, resample=resample, require=require)
    df1.to_csv(csv_path, header=True)


def to_netcdf(self, nc_path):
    """
    Function to export a copy of a hydro class object to a netcdf file.
    """
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
        crs1 = geo1.crs
        if isinstance(crs1, str):
            crs2 = convert_crs(crs1, 'proj4_dict')
        elif isinstance(crs1, dict):
            crs2 = crs1.copy()
        crs2.update({i: str(crs2[i]) for i in crs2 if type(crs2[i]) == bool})
        ds0['geo_loc_x'].attrs = crs2
        ds0['geo_loc_y'].attrs = crs2


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

    Parameters
    ----------
    shp_path : str
        The shapefile path.

    Returns
    -------
    None
    """

    if not hasattr(self, 'geo_loc'):
        raise ValueError('Object has no geo locations!')

    ### Prepare output
    geo1 = self.geo_loc.copy()
    sites = self.sites
    mtypes_sites = self.mtypes_sites.copy()

    if len(mtypes_sites.keys()) > 1:
        t1 = {i: in1d(sites, list(mtypes_sites[i])).astype(str).tolist() for i in mtypes_sites}
        df1 = DataFrame(t1, index=sites)
        geo2 = concat([geo1, df1], axis=1).reset_index()
        geo2.to_file(shp_path)
    else:
        self._base_stats_fun()
        stats = self._base_stats.copy()
        stats.index = stats.index.droplevel('mtype')
        stats['start_time'] = stats['start_time'].astype(str)
        stats['end_time'] = stats['end_time'].astype(str)
        geo2 = concat([geo1, stats], axis=1).reset_index()
        geo2.to_file(shp_path)













