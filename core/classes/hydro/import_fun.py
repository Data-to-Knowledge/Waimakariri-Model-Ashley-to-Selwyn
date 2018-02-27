# -*- coding: utf-8 -*-
"""
Functions for the hydro class for importing data.
"""

from pandas import DataFrame, Series, DatetimeIndex, to_datetime, MultiIndex, concat, to_numeric, read_csv
from numpy import in1d, where
from geopandas import GeoDataFrame, GeoSeries, read_file
from collections import defaultdict
from core.classes.hydro.base import all_mtypes
from xarray import open_dataset
from shapely.wkt import loads
from shapely.geometry import Point
from core.ecan_io.mssql import rd_sql_ts, rd_sql
from core.spatial.vector import sel_sites_poly


########################################################
#### Time series data

### Primary import functions
#h1.add_data(df, time='time', sites='site', mtypes='mtype', values='data', dformat='long')

def add_data(self, data, time=None, sites=None, mtypes=None, values=None, dformat=None, add=True):
    """
    The general function to add time series data to the hydro class object.

    Input data can be either in 'wide' or 'long' dformat.

    The 'wide' dformat is where the columns are the sites or both mtypes and sites as a MultiIndex column. 'time' should either be None if data.index is a DateTimeIndex or a DateTimeIndex. 'sites' should be None. If data has a MultiIndex column of mtypes and sites, they must be arranged in that order (top: mtypes, bottom: sites).

    The 'long' dformat should be a DataFrame with four columns: One column with a DatetimeIndex for 'time', one column with the site values as 'sites', one column with the mtypes values as 'mtypes', and one column with the data values as 'values'.

    Parameters
    ----------
    data : DataFrame or type that can be coerced to a DataFrame
        A Pandas DataFrame in either long or wide format or an object type that can be coerced to a DataFrame via DataFrame(data) (see Pandas).
    time : DateTimeIndex, str, or None
        The time index reference. Depends on dformat.
    sites : str, int, or None
        The sites reference. Depends on dformat.
    mtypes : str, int, or None
        The mtypes reference. Depends on dformat.
    values : str or None
        Only needed for dformat = 'long'. The column name for the data values.
    dformat : 'wide' or 'long'
        The format of the data table to import.
    add : bool
        Should new data be appended to the existing object? If False, returns a new object.

    Returns
    -------
    HydroPandas
    """

    ### Convert input data to long format for consumption

    if dformat is None:
        raise ValueError("dformat must be specified and must be either 'long' or 'wide'.")
    if not isinstance(data, DataFrame):
        data = DataFrame(data)
    if data.empty:
        raise ValueError('DataFrame is empty, no data was passed')
    if dformat is 'wide':
        if isinstance(data.columns, MultiIndex):
            if not isinstance(data.index, DatetimeIndex):
                raise ValueError("A MultiIndex column DataFrame must have a DateTimeIndex as the row index.")
            d1 = data.copy()
            d1.columns.names = ['mtype', 'site']
            d1.columns.set_levels(to_numeric(d1.columns.levels[1], errors='ignore', downcast='integer'), level='site', inplace=True)
            d1.index.name = 'time'
            d2 = d1.stack(level=['mtype', 'site'])
            input1 = d2.reorder_levels(['mtype', 'site', 'time'])
        else:
            if isinstance(data.index, DatetimeIndex):
                    d1 = data.copy()
                    d1.index.name = 'time'
            else:
                if isinstance(time, (str, int)):
                    data.loc[:, time] = to_datetime(data.loc[:, time])
                    d1 = data.set_index(time)
                else:
                    d1 = data.set_index(to_datetime(time))
                d1.index.name = 'time'
            if isinstance(mtypes, (str, int)):
                d1.columns.name = 'site'
                d1.columns = to_numeric(d1.columns, errors='ignore', downcast='integer')
                d1['mtype'] = mtypes
                d2 = d1.set_index('mtype', append=True)
                input1 = d2.stack().reorder_levels(['mtype', 'site', 'time'])
    elif dformat is 'long':
        if not isinstance(mtypes, (str, int)) or not isinstance(sites, (str, int)) or not isinstance(time, (str, int)) or not isinstance(values, (str, int)):
            raise ValueError("Long format DataFrame must have 4 columns each corresponding to 'mtypes', 'sites', 'time', and 'values'.")
        else:
            d1 = data.copy()
            cols = d1.columns
            if len(cols) < 4:
                set1 = [mtypes, sites]
                w1 = where(Series(set1).isin(cols))[0][0]
                if w1 == 0:
                    d1['site'] = sites
                    sites = 'site'
                else:
                    d1['mtype'] = mtypes
                    mtypes = 'mtype'

            if not isinstance(d1.loc[:, time], DatetimeIndex):
                d1.loc[:, time] = to_datetime(d1.loc[:, time])
            d1.rename(columns={mtypes: 'mtype', sites: 'site', time: 'time'}, inplace=True)
            input1 = d1.set_index(['mtype', 'site', 'time'])[values]
            input1.index.set_levels(to_numeric(input1.index.levels[1], errors='ignore'), level='site', inplace=True)

    ### Run additional checks and append to existing data
    ## Remove data with NaN
    if any(input1.isnull()):
        input2 = input1[input1.notnull()].reset_index()
        input1 = input2.set_index(['mtype', 'site', 'time'])[values]

    ## Check mtypes
    new_mtypes = input1.index.levels[0]
    mtypes_bool = in1d(new_mtypes, all_mtypes)
    if any(~mtypes_bool):
        sel_mtypes = new_mtypes[mtypes_bool]
        input1 = input1.loc(axis=0)[sel_mtypes, :, :]

    input1.name = 'data'

    ## Remove base_stats if it exists
    if hasattr(self, '_base_stats'):
        delattr(self, '_base_stats')

    ## Add data to hydro object
    if add:
        if hasattr(self, 'data'):
            setattr(self, 'data', input1.combine_first(self.data).sort_index())
        else:
            setattr(self, 'data', input1.sort_index())

        ## Reassign attributes
        _import_attr(self)

        ## Return
        return(self)
    else:
        new1 = self.copy()
        setattr(new1, 'data', input1.sort_index())
        _import_attr(new1)
        if hasattr(new1, 'site_geo_attr'):
            new1.get_site_geo_attr()
        if hasattr(new1, 'geo_loc'):
            new1.add_geo_loc(new1.geo_loc, check=False)
        if hasattr(new1, 'geo_catch'):
            new1.add_geo_catch(new1.geo_catch, check=False)
        return(new1)


### Import helper functions


### Functions applied after time series import to assign new attributes

def _import_attr(self):
    ## Assign the list of available mtypes
    self.mtypes = self.data.index.levels[0].tolist()

    ## Assign the list of available sites
    self.sites = self.data.index.levels[1].tolist()

    ## Assign a dict of sites to mtypes
    mtypes1 = self.data.index.get_level_values('mtype')
    sites1 = self.data.index.get_level_values('site')
    zipped = zip(mtypes1, sites1)
    ms = defaultdict(set)
    for a, b in zipped:
        ms[a].add(b)
    setattr(self, 'mtypes_sites', dict(ms))

    sm = defaultdict(set)
    for a, b in zipped:
        sm[b].add(a)
    setattr(self, 'sites_mtypes', dict(sm))

    ## Assign a blank site attribute table
#    if hasattr(self, 'site_attr'):
#        sites_needed = set(self.sites).difference(self.site_attr.index.tolist())
#        new_sites = DataFrame(nan, index=sites_needed, columns=site_attr_lst)
#        combine = concat([self.site_attr, new_sites])
#        combine.index.name = 'site'
#        setattr(self, 'site_attr', combine)
#    else:
#        site_attr = DataFrame(nan, index=self.sites, columns=site_attr_lst)
#        site_attr.index.name = 'site'
#        setattr(self, 'site_attr', site_attr)

###################################################
#### Combine two hydro class objects together

def combine(self, hp):
    """
    Function to combine two hydro class objects.
    """
    new1 = self.add_data(hp.data.reset_index(), 'time', 'site', 'mtype', 'data', dformat='long')
    if hasattr(hp, 'site_attr'):
        new1.add_site_attr(hp.site_attr)
    if hasattr(hp, 'geo_loc'):
        new1.add_geo_loc(hp.geo_loc, check=False)
    if hasattr(hp, 'geo_catch'):
        new1.add_geo_catch(hp.geo_catch, check=False)
    return(new1)


####################################################
#### Geometry data

### Import geometry with sites

def _add_geo_data(self, geo, obj_name):
    """
    Function to read in a GeoDataFrame to the hydro class. The GeoDataFrame must only contain a geometry column and the index must be the site names. Or it can be a shapefile where the first column has the site names.
    """
    if hasattr(self, obj_name):
        pre_geo = getattr(self, obj_name)
        geo_sites = pre_geo.index.tolist()
    else:
        geo_sites = []
    sites = self.sites
    geo_sites_check = ~Series(sites).isin(geo_sites)
    if (len(sites) < len(geo_sites)) & all(~geo_sites_check):
        pre_geo = pre_geo.loc[sites]
        setattr(self, obj_name, pre_geo)
    if any(geo_sites_check):
        sites_sel = Series(sites)[geo_sites_check]
        if isinstance(geo, (GeoDataFrame, str)):
            if isinstance(geo, str):
                geo0 = read_file(geo)
                geo1 = geo0.set_index(geo0.columns[0])
            else:
                geo1 = geo
            site_sel_check = Series(sites_sel).isin(geo1.index)
            sites_sel2 = sites_sel[site_sel_check]
            if any(site_sel_check):
                if hasattr(self, obj_name):
                    setattr(self, obj_name, concat([pre_geo, geo1.loc[sites_sel2, :]]))
                else:
                    setattr(self, obj_name, geo1.loc[sites_sel2, :])
                getattr(self, obj_name).index.name = 'site'
            else:
                raise ValueError('Input geo data must have locations for some of the mtype sites.')
        else:
            raise ValueError('Input data must be a GeoDataFrame or a string path to a shapefile!')


def add_geo_loc(self, geo, check=True):
    _add_geo_data(self, geo, 'geo_loc')
    if check:
        _check_geo_sites(self, 'geo_loc')


def add_geo_catch(self, geo, check=True):
    if hasattr(self, 'geo_loc'):
        from core.classes.hydro.indexing import _comp_by_catch
        _add_geo_data(self, geo, 'geo_catch')
        _comp_by_catch(self)
        if check:
            _check_crs(self)
            _check_geo_sites(self, 'geo_catch')
    else:
        raise ValueError('Add the site locations first before adding the catchment polygons.')

### Checks

## Identify missing site locations
def missing_geo_loc_sites(self):
    if hasattr(self, 'geo_loc'):
        geo_sites = self.geo_loc.index.tolist()
        sites = self.sites
        mis_sites = Series(sites)[~Series(sites).isin(geo_sites)].tolist()
        self._mis_geo_sites = mis_sites
        return(mis_sites)
    else:
        raise ValueError('You are missing all of the sites!')

## Make sure that the geo objects have the same crs (reproject if needed)
def _check_crs(self):
    geo_crs = {i: getattr(self, i).crs for i in self.__dict__.keys() if 'geo' in i}
    if len(geo_crs) > 1:
        geo_loc_crs = getattr(self, 'geo_loc').crs
        for i in geo_crs.keys():
            if not i == 'geo_loc':
                geo1 = getattr(self, i)
                geo2 = geo1.to_crs(crs=geo_loc_crs)
                setattr(self, i, geo2)

## Check if all sites have been found for the geo data
def _check_geo_sites(self, geo_name):
    all_sites = set(self.sites)
    geo_sites = getattr(self, geo_name).index.tolist()
    mis1 = all_sites.difference(geo_sites)
    if mis1:
        print('Missing sites ' + str(list(mis1)))
    else:
        print('Found all of the sites!')

####################################################
#### Populate site attributes

def add_site_attr(self, df):
    """
    Function to add in site attributes. df must be a dataframe with the sites as an index.
    """
    if isinstance(df, DataFrame):
        sites = self.sites
        if all(Series(sites).isin(df.index)):
            new_attr1 = df.loc[df.index.isin(sites), site_attr_lst]
        elif any(Series(sites).isin(df.index)):
            new_attr = df.loc[df.index.isin(sites), site_attr_lst]
            new_attr1 = new_attr.combine_first(self.site_attr)
        else:
            raise ValueError('No new sites!')
        setattr(self, 'site_attr', new_attr1)


####################################################
#### File reading functions (csv and netcdf)


def rd_csv(self, csv_path, time=None, sites=None, mtypes=None, values=None, dformat=None, multiindex=False, skiprows=0):
    """
    Simple function to read in time series data and make it regular if needed.
    """

    if multiindex:
        ts = read_csv(csv_path, parse_dates=True, infer_datetime_format=True, dayfirst=True, skiprows=skiprows, header=[0, 1], index_col=0)
        ts.columns.names = ['mtype', 'site']
        ts.index.name = 'time'
    else:
        ts = read_csv(csv_path, parse_dates=[time], infer_datetime_format=True, dayfirst=True, skiprows=skiprows, header='infer', index_col=0)
    self.add_data(ts, time=time, sites=sites, mtypes=mtypes, values=values, dformat=dformat)
    return(self)


def rd_netcdf(self, nc_path):
    """
    Function to read a netcdf file (.nc) that was an export from a hydro class.
    """

    with open_dataset(nc_path) as ds1:

        ### Load in the geometry data
        if any(in1d(ds1.data_vars.keys(), 'geo_catch_wkt')):
    #        geo_catch_cols = [i for i in ds1.keys() if 'geo_catch' in i]
            df_catch = ds1['geo_catch_wkt'].to_dataframe()
            df_catch.columns = df_catch.columns.str.replace('geo_catch_', '')
            geo1 = [loads(x) for x in df_catch.wkt]
            gdf_catch = GeoDataFrame(df_catch.drop('wkt', axis=1), geometry=geo1, crs=dict(ds1.attrs['crs']))
            gdf_catch.index.name = 'site'
            ds1 = ds1.drop('geo_catch_wkt')

        if any(in1d(ds1.data_vars.keys(), 'geo_loc_x')):
            geo_loc_cols = [i for i in ds1.keys() if 'geo_loc' in i]
            df_loc = ds1[geo_loc_cols].to_dataframe()
            df_loc.columns = df_loc.columns.str.replace('geo_loc_', '')
            geo2 = [Point(xy) for xy in zip(df_loc.x, df_loc.y)]
            crs1 = dict(ds1['geo_loc_x'].attrs.copy())
            crs1.update({i: bool(crs1[i]) for i in crs1 if crs1[i] in ['True', 'False']})
            gdf_loc = GeoDataFrame(df_loc.drop(['x', 'y'], axis=1), geometry=geo2, crs=crs1)
            gdf_loc.index.name = 'site'
            ds1 = ds1.drop(geo_loc_cols)

        ### Load in the site attribute data
        if any(in1d(ds1.data_vars.keys(), 'site_attr')):
            site_attr_cols = [i for i in ds1.keys() if 'site_attr' in i]
            site_attr1 = ds1[site_attr_cols].to_dataframe()
            site_attr1.columns = site_attr1.columns.str.replace('site_attr_', '')
            site_attr1.index.name = 'site'
            ds1 = ds1.drop(site_attr_cols)

        ### Load in the ts data
        df1 = ds1[['site', 'time', 'data']].to_dataframe().reset_index()

        self.add_data(df1, 'time', 'site', 'mtype', 'data', 'long')

        ### Add in the earlier attributes
        if 'site_attr1' in locals():
            self.add_site_attr(site_attr1)
        if 'gdf_loc' in locals():
            self.add_geo_loc(gdf_loc, check=False)
        if 'gdf_catch' in locals():
            self.add_geo_catch(gdf_catch, check=False)

        ### return
    #    ds1.close()
        return(self)


##############################################
#### Functions to read from databases

### mssql


def _rd_hydro_mssql(self, server, database, table, mtype, date_col, site_col, data_col, qual_col=None, sites=None, from_date=None, to_date=None, qual_codes=None, add_where=None, min_count=None, resample_code=None, period=1, fun='mean', val_round=3):
    """
    Function to import data from a MSSQL database. Specific columns can be selected and specific queries within columns can be selected. Requires the pymssql package.

    Parameters
    ----------
    server : str
        The server name. e.g.: 'SQL2012PROD03'
    database : str
        The specific database within the server. e.g.: 'LowFlows'
    table : str
        The specific table within the database. e.g.: 'LowFlowSiteRestrictionDaily'
    mtype : str
        The measurement type according to the Hydro mtypes.
    date_col : str
        The DateTime column name.
    site_col : str
        The site name column.
    data_col : str
        The data column.
    qual_col : str
        The quality code column.
    sites : list
        List of sites.
    from_date : str
        From date for data in the format '2000-01-01'.
    to_date : str
        To date for the data in the above format.
    qual_codes : list of int
        List of quality codes.
    add_where : str
        An additional SQL query to be added.
    min_count : int
        The minimum number of data values per site for data extraction.
    resample_code : str or None
        The Pandas time series resampling code. e.g. 'D' for day, 'W' for week, 'M' for month, etc.
    period : int
        The number of resampling periods. e.g. period = 2 and resample = 'D' would be to resample the values over a 2 day period.
    fun : str
        The resampling function. i.e. mean, sum, count, min, or max. No median yet...
    val_round : int
        The number of decimals to round the values.

    Return
    ------
    Hydro
    """

    ## Prepare sql_dict for sites and qual_cdes
    site_qual_dict = {}

    if isinstance(qual_codes, list):
        site_qual_dict.update({qual_col: qual_codes})

    if isinstance(sites, list):
        sites = [str(i) for i in sites]
        site_qual_dict.update({site_col: sites})

    sql_dict = {'server': server, 'database': database, 'table': table, 'groupby_cols': site_col, 'date_col': date_col, 'values_cols': data_col, 'resample_code': resample_code, 'period': period, 'fun': fun, 'val_round': val_round, 'where_col': site_qual_dict, 'from_date': from_date, 'to_date': to_date, 'min_count': min_count}

    df = rd_sql_ts(**sql_dict).reset_index()

    ## Rename columns
    df.columns = ['site', 'time', 'data']

    ## Remove spaces in site names and duplicate data
    df.loc[:, 'site'] = df.loc[:, 'site'].astype(str).str.replace(' ', '').str.upper()
    df = df.drop_duplicates(['site', 'time'])

    df['mtype'] = mtype

    self.add_data(df, 'time', 'site', 'mtype', 'data', 'long')
    return(self)


def _rd_hydro_geo_mssql(self, server, database, table, geo_dict):
    """
    Function to select sites based on the geo attributes.
    """

    sites1 = rd_sql(server, database, table, 'site', geo_dict)
    sites2 = sites1.site.astype(str).values.tolist()
    return(sites2)


def _proc_hydro_sql(self, sites_sql_fun, mtype_dict, mtype, sites=None, from_date=None, to_date=None, qual_codes=None, min_count=None, buffer_dis=0, resample_code=None, period=1, fun='mean'):
    """
    Convenience function for reading in mssql data from standardized hydro tables.
    """

    if isinstance(sites, GeoDataFrame):
        loc1 = sites_sql_fun()
        sites1 = sel_sites_poly(loc1, sites, buffer_dis).index.astype(str)
    else:
        sites1 = Series(sites).astype(str)

    h1 = self.copy()
    if isinstance(mtype_dict, (list, tuple)):
        for i in range(len(mtype_dict)):
            site1 = mtype_dict[i]['site_col']

            sites_stmt = 'select distinct ' + site1 + ' from ' + mtype_dict[i]['table']
            sites2 = rd_sql(mtype_dict[i]['server'], mtype_dict[i]['database'], stmt=sites_stmt).astype(str)[site1]
            sites3 = sites2[sites2.isin(sites1)].astype(str).tolist()
            if not sites3:
                raise ValueError('No sites in database')
            if mtype_dict[i]['qual_col'] is None:
                qual_codes = None
            h1 = h1._rd_hydro_mssql(sites=sites3, mtype=mtype, from_date=from_date, to_date=to_date, qual_codes=qual_codes, min_count=min_count, resample_code=resample_code, period=period, fun=fun, **mtype_dict[i])
    elif isinstance(mtype_dict, dict):
        site1 = mtype_dict['site_col']

        sites_stmt = 'select distinct ' + site1 + ' from ' + mtype_dict['table']
        sites2 = rd_sql(mtype_dict['server'], mtype_dict['database'], stmt=sites_stmt).astype(str)[site1]
        sites3 = sites2[sites2.isin(sites1)].astype(str).tolist()
        if not sites3:
                raise ValueError('No sites in database')
        if mtype_dict['qual_col'] is None:
            qual_codes = None
        h1 = h1._rd_hydro_mssql(sites=sites3, mtype=mtype, from_date=from_date, to_date=to_date, qual_codes=qual_codes, min_count=min_count, resample_code=resample_code, period=period, fun=fun, **mtype_dict)
    elif callable(mtype_dict):
        h1 = mtype_dict(h1, sites=sites1, mtype=mtype, from_date=from_date, to_date=to_date, min_count=min_count)

    return(h1)







