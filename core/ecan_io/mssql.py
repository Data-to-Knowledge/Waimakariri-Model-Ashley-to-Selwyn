# -*- coding: utf-8 -*-
"""
Functions for importing mssql data.
"""
from geopandas import GeoDataFrame
from core.misc.misc import save_df
from pandas import to_numeric, to_datetime, Timestamp, read_sql, Series, concat, merge, DataFrame
from shapely.wkt import loads
from pycrs.parser import from_epsg_code
from pymssql import connect
import traceback
from datetime import datetime
#from core.misc.misc import select_sites
#from core.spatial.vector import xy_to_gpd, sel_sites_poly
#from numpy import nan, ceil


def rd_sql(server, database, table=None, col_names=None, where_col=None, where_val=None, where_op='AND', geo_col=False,
           from_date=None, to_date=None, date_col=None, rename_cols=None, stmt=None, export_path=None):
    """
    Function to import data from an MSSQL database. Requires the pymssql package.

    Parameters
    ----------
    server : str
        The server name. e.g.: 'SQL2012PROD03'
    database : str
        The specific database within the server. e.g.: 'LowFlows'
    table : str
        The specific table within the database. e.g.: 'LowFlowSiteRestrictionDaily'
    col_names : list of str
        The column names that should be retrieved. e.g.: ['SiteID', 'BandNo', 'RecordNo']
    where_col : str or dict
        Must be either a string with an associated where_val list or a dictionary of strings to lists.'. e.g.: 'SnapshotType' or {'SnapshotType': ['value1', 'value2']}
    where_val : list
        The WHERE query values for the where_col. e.g. ['value1', 'value2']
    where_op : str
        If where_col is a dictionary and there are more than one key, then the operator that connects the where statements must be either 'AND' or 'OR'.
    geo_col : bool
        Is there a geometry column in the table?.
    from_date : str
        The start date in the form '2010-01-01'.
    to_date : str
        The end date in the form '2010-01-01'.
    date_col : str
        The SQL table column that contains the dates.
    stmt : str
        Custom SQL statement to be directly passed to the database table. This will ignore all prior arguments except server and database.
    export_path : str
        The export path for a csv file if desired. If None, then nothing is exported.

    Returns
    -------
    DataFrame
    """

    ## Create where statements
    if stmt is None:

        if table is None:
            raise ValueError('Must at least provide input for server, database, and table.')

        if col_names is not None:
            if isinstance(col_names, str):
                col_names = [col_names]
            col_names1 = ['[' + i.encode('ascii', 'ignore') + ']' for i in col_names]
            col_stmt = ', '.join(col_names1)
        else:
            col_stmt = '*'

        where_lst = sql_where_stmts(where_col=where_col, where_val=where_val, where_op=where_op, from_date=from_date,
                                    to_date=to_date, date_col=date_col)

        if isinstance(where_lst, list):
            stmt1 = "SELECT " + col_stmt + " FROM " + table + " where " + " and ".join(where_lst)
        else:
            stmt1 = "SELECT " + col_stmt + " FROM " + table

    elif isinstance(stmt, str):
        stmt1 = stmt

    else:
        raise ValueError('stmt must either be an SQL string or None.')

    ## Create connection to database and execute sql statement
    conn = connect(server, database=database)
    df = read_sql(stmt1, conn)
    conn.close()
    if rename_cols is not None:
        df.columns = rename_cols

    ## Read in geometry if required
    if geo_col & (stmt is not None):
        geometry, proj = rd_sql_geo(server=server, database=database, table=table, where_lst=where_lst)
        df = GeoDataFrame(df, geometry=geometry, crs=proj)

    ## save and return
    if export_path is not None:
        save_df(df, export_path, index=False)

    return (df)


def rd_sql_ts(server, database, table, groupby_cols, date_col, values_cols, resample_code=None, period=1, fun='mean',
              val_round=3, where_col=None, where_val=None, where_op='AND', from_date=None, to_date=None, min_count=None,
              export_path=None):
    """
    Function to specifically read and possibly aggregate time series data stored in MSSQL tables. Returns a MultiIndex DataFrame.

    Parameters
    ----------
    server : str
        The server name. e.g.: 'SQL2012PROD03'
    database : str
        The specific database within the server. e.g.: 'LowFlows'
    table : str
        The specific table within the database. e.g.: 'LowFlowSiteRestrictionDaily'
    groupby_cols : str or list of str
        The columns in the SQL table to grouped and returned with the time series data.
    date_col : str
        The date column in the table.
    values_cols : str or list of str
        The column(s) of the value(s) that should be resampled.
    resample_code : str or None
        The Pandas time series resampling code. e.g. 'D' for day, 'W' for week, 'M' for month, etc.
    period : int
        The number of resampling periods. e.g. period = 2 and resample = 'D' would be to resample the values over a 2 day period.
    fun : str
        The resampling function. i.e. mean, sum, count, min, or max. No median yet...
    val_round : int
        The number of decimals to round the values.
    where_col : str or dict
        Must be either a string with an associated where_val list or a dictionary of strings to lists.'. e.g.: 'SnapshotType' or {'SnapshotType': ['value1', 'value2']}
    where_val : list
        The WHERE query values for the where_col. e.g. ['value1', 'value2']
    where_op : str
        If where_col is a dictionary and there are more than one key, then the operator that connects the where statements must be either 'AND' or 'OR'.
    from_date : str
        The start date in the form '2010-01-01'.
    to_date : str
        The end date in the form '2010-01-01'.
    min_count : int
        The minimum number of values required to return groupby_cols. Only works when groupby_cols and vlue_cols are str.
    export_path : str
        The export path for a csv file if desired. If None, then nothing is exported.

    Returns
    -------
    DataFrame
        Pandas DataFrame with MultiIndex of groupby_cols and date_col
    """

    ## Create where statement
    where_lst = sql_where_stmts(where_col=where_col, where_val=where_val, where_op=where_op, from_date=from_date,
                                to_date=to_date, date_col=date_col)

    ## Create ts statement and append earlier where statement
    if isinstance(groupby_cols, str):
        groupby_cols = [groupby_cols]
    col_names1 = ['[' + i.encode('ascii', 'ignore') + ']' for i in groupby_cols]
    col_stmt = ', '.join(col_names1)

    ## Create sql stmt
    sql_stmt1 = sql_ts_agg_stmt(table, groupby_cols=groupby_cols, date_col=date_col, values_cols=values_cols,
                                resample_code=resample_code, period=period, fun=fun, val_round=val_round,
                                where_lst=where_lst)

    ## Create connection to database
    conn = connect(server, database=database)

    ## Make minimum count selection
    if (min_count is not None) & isinstance(min_count, int) & (len(groupby_cols) == 1):
        cols_count_str = ', '.join([groupby_cols[0], 'count(' + values_cols + ') as count'])
        if isinstance(where_lst, list):
            stmt1 = "SELECT " + cols_count_str + " FROM " + "(" + sql_stmt1 + ") as agg" + " GROUP BY " + col_stmt + " HAVING count(" + values_cols + ") >= " + str(
                min_count)
        else:
            stmt1 = "SELECT " + cols_count_str + " FROM " + table + " GROUP BY " + col_stmt + " HAVING count(" + values_cols + ") >= " + str(
                min_count)

        up_sites = read_sql(stmt1, conn)[groupby_cols[0]].tolist()
        up_sites = [str(i) for i in up_sites]

        if not up_sites:
            raise ValueError('min_count filtered out all sites.')

        if isinstance(where_col, str):
            where_col = {where_col: where_val}

        where_col.update({groupby_cols[0]: up_sites})
        where_lst = sql_where_stmts(where_col, where_op=where_op, from_date=from_date, to_date=to_date,
                                    date_col=date_col)

        ## Create sql stmt
        sql_stmt1 = sql_ts_agg_stmt(table, groupby_cols=groupby_cols, date_col=date_col, values_cols=values_cols,
                                    resample_code=resample_code, period=period, fun=fun, val_round=val_round,
                                    where_lst=where_lst)

    ## Create connection to database and execute sql statement
    df = read_sql(sql_stmt1, conn)
    conn.close()

    ## Check to see if any data was found
    if df.empty:
        raise ValueError('No data was found in the database for the parameters given.')

    ## set the index
    groupby_cols.append(date_col)
    df1 = df.set_index(groupby_cols).sort_index()

    ## Save and return
    if export_path is not None:
        save_df(df1, export_path)

    return (df1)


#def rd_sw_rain_geo(sites=None):
#    if sites is not None:
#        site_geo = rd_sql('SQL2012PROD05', 'Bgauging', 'RSITES',
#                          col_names=['SiteNumber', 'River', 'SiteName', 'NZTMX', 'NZTMY'], where_col='SiteNumber',
#                          where_val=sites)
#    else:
#        site_geo = rd_sql('SQL2012PROD05', 'Bgauging', 'RSITES',
#                          col_names=['SiteNumber', 'River', 'SiteName', 'NZTMX', 'NZTMY'])
#
#    site_geo.columns = ['site', 'river', 'name', 'NZTMX', 'NZTMY']
#    site_geo.loc[:, 'site'] = to_numeric(site_geo.loc[:, 'site'], errors='ignore')
#
#    site_geo2 = xy_to_gpd(df=site_geo, id_col=['site', 'river', 'name'], x_col='NZTMX', y_col='NZTMY')
#    site_geo3 = site_geo2.loc[site_geo2.site > 0, :]
#    site_geo3.loc[:, 'site'] = site_geo3.loc[:, 'site'].astype('int32')
#    site_geo3['river'] = site_geo3.river.apply(lambda x: x.title())
#    site_geo3['name'] = site_geo3.name.apply(lambda x: x.title())
#    site_geo3['name'] = site_geo3.name.apply(lambda x: x.replace(' (Recorder)', ''))
#    site_geo3['name'] = site_geo3.name.apply(lambda x: x.replace('Sh', 'SH'))
#    site_geo3['name'] = site_geo3.name.apply(lambda x: x.replace('Ecs', 'ECS'))
#    return (site_geo3.set_index('site'))


def write_sql(df, server, database, table, dtype_dict, primary_keys=None, foreign_keys=None, foreign_table=None, create_table=True, drop_table=False, del_rows_dict=None, output_stmt=False):
    """
    Function to write pandas dataframes to mssql server tables. Must have write permissions to database!

    Parameters
    ----------
    df : DataFrame
        DataFrame to be saved.
    server : str
        The server name. e.g.: 'SQL2012PROD03'
    database : str
        The specific database within the server. e.g.: 'LowFlows'
    table : str
        The specific table within the database. e.g.: 'LowFlowSiteRestrictionDaily'
    dtype_dict : dict of str
        Dictionary of df columns to the associated sql data type. Examples below.
    primary_keys : str or list of str
        Index columns to define uniqueness in the data structure.
    foreign_keys : str or list of str
        Columns to link to another table in the same database.
    foreign_table: str
        The table in the same database with the identical foreign key(s).
    create_table : bool
        Should a new table be created or should it be appended to an existing table?
    drop_table : bool
        If the table already exists, should it be dropped?
    output_stmt : bool
        Should the SQL statements be outputted to a dictionary?

    Returns
    -------
    if output_stmt is True then dict else None

    dtype strings for matching python data types to SQL
    ---------------------------------------------------
    str: 'VARCHAR(19)'

    date: 'DATE'

    datetime: "DATETIME'

    float: 'NUMERIC(10, 1)' or 'FLOAT'

    int: 'INT'
    """

    #### Parameters and functions
    py_sql = {'NUMERIC': float, 'DATE': str, 'DATETIME': str, 'INT': 'int32', 'smallint': 'int32', 'date': str, 'varchar': str, 'float': float, 'datetime': str, 'VARCHAR': str, 'FLOAT': float, 'smalldatetime': str, 'decimal': float}

    def chunker(seq, size):
        return ([seq[pos:pos + size] for pos in range(0, len(seq), size)])

    #### Make sure the df has the correct dtypes
    if len(dtype_dict) != len(df.columns):
        raise ValueError('dtype_dict must have the same number of keys as columns in df.')
    if not all(df.columns.isin(dtype_dict.keys())):
        raise ValueError('dtype_dict must have the same column names as the columns in the df.')

    df1 = df.copy()

    for i in df.columns:
        dtype1 = dtype_dict[i]
        if (dtype1 == 'DATE') | (dtype1 == 'DATETIME'):
            time1 = to_datetime(df[i]).dt.isoformat(' ')
            df1.loc[:, i] = time1
        elif 'VARCHAR' in dtype1:
            try:
                df1.loc[:, i] = df.loc[:, i].astype(str).str.replace('\'', '')
            except:
                df1.loc[:, i] = df.loc[:, i].str.encode('utf-8', 'ignore').str.replace('\'', '')
        elif 'NUMERIC' in dtype1:
            df1.loc[:, i] = df.loc[:, i].astype(float)
        elif 'decimal' in dtype1:
            df1.loc[:, i] = df.loc[:, i].astype(float)
        elif not dtype1 in py_sql.keys():
            raise ValueError('dtype must be one of ' + str(py_sql.keys()))
        else:
            df1.loc[:, i] = df.loc[:, i].astype(py_sql[dtype1])

    #### Convert df to set of tuples to be ingested by sql
    list1 = df1.values.tolist()
    tup1 = [str(tuple(i)) for i in list1]
    tup2 = chunker(tup1, 1000)

    #### Primary keys
    if isinstance(primary_keys, str):
        primary_keys = [primary_keys]
    if isinstance(primary_keys, list):
        key_stmt = ", Primary key (" + ", ".join(primary_keys) + ")"
    else:
        key_stmt = ""

    #### Foreign keys
    if isinstance(foreign_keys, str):
        foreign_keys = [foreign_keys]
    if isinstance(foreign_keys, list):
        fkey_stmt = ", Foreign key (" + ", ".join(foreign_keys) + ") " + "References " + foreign_table + "(" + ", ".join(foreign_keys) + ")"
    else:
        fkey_stmt = ""

    #### Initial create table and insert statements
    d1 = [str(i) + ' ' + dtype_dict[i] for i in df.columns]
    d2 = ', '.join(d1)
    tab_create_stmt = "create table " + table + " (" + d2 + key_stmt + fkey_stmt + ")"
    columns1 = str(tuple(df1.columns.tolist())).replace('\'', '')
    insert_stmt1 = "insert into " + table + " " + columns1 + " values "

    conn = connect(server, database=database)
    cursor = conn.cursor()
    stmt_dict = {}
    try:

        #### Drop table if it exists
        if drop_table:
            drop_stmt = "IF OBJECT_ID(" + str([str(table)])[1:-1] + ", 'U') IS NOT NULL DROP TABLE " + table
            cursor.execute(drop_stmt)
            conn.commit()
            stmt_dict.update({'drop_stmt': drop_stmt})

        #### Create table in database
        if create_table:
            cursor.execute(tab_create_stmt)
            conn.commit()
            stmt_dict.update({'tab_create_stmt': tab_create_stmt})

        #### Delete rows
        if isinstance(del_rows_dict, dict):
            del_where_list = sql_where_stmts(**del_rows_dict)
            del_rows_stmt = "DELETE FROM " + table + " WHERE " + " AND ".join(del_where_list)
            cursor.execute(del_rows_stmt)
            conn.commit()

        #### Insert data into table
        for i in range(len(tup2)):
            rows = ",".join(tup2[i])
            insert_stmt2 = insert_stmt1 + rows
            cursor.execute(insert_stmt2)
            stmt_dict.update({'insert' + str(i+1): insert_stmt2})
        conn.commit()

        #### Close everything!
        cursor.close()
        conn.close()

        if output_stmt:
            return stmt_dict
    except Exception as err:
        conn.rollback()
        conn.close()
        raise err


def sql_where_stmts(where_col=None, where_val=None, where_op='AND', from_date=None, to_date=None, date_col=None):
    """
    Function to take various input parameters and convert them to a list of where statements for SQL.

    Parameters
    ----------
    where_col : str or dict
        Either a str with an associated where_val list or a dictionary of string keys to list values. If a str, it should represent the table column associated with the 'where' condition.
    where_val : list or None
        If where_col is a str, then where_val must be a list of associated condition values.
    where_op : str of either 'AND' or 'OR'
        The binding operator for the where conditions.
    from_date : str or None
        The start date in the form '2010-01-01'.
    to_date : str or None
        The end date in the form '2010-01-01'.
    date_col : str or None
        The SQL table column that contains the dates.

    Returns
    -------
    list of str or None
        Returns a list of str where conditions to be passed to an SQL execution function. The function needs to bind it with " where " + " and ".join(where_lst)
    """

    if where_col is not None:
        if isinstance(where_col, str) & isinstance(where_val, list):
            #            if len(where_val) > 10000:
            #                raise ValueError('The number of values in where_val cannot be over 10000 (or so). MSSQL limitation. Break them into smaller chunks.')
            where_val = [str(i) for i in where_val]
            where_stmt = [str(where_col) + ' IN (' + str(where_val)[1:-1] + ')']
        elif isinstance(where_col, dict):
            where_stmt = [i + " IN (" + str([str(j) for j in where_col[i]])[1:-1] + ")" for i in where_col]
        else:
            raise ValueError(
                'where_col must be either a string with an associated where_val list or a dictionary of string keys to list values.')
    else:
        where_stmt = []

    if isinstance(from_date, str):
        from_date1 = to_datetime(from_date, errors='coerce')
        if isinstance(from_date1, Timestamp):
            from_date2 = from_date1.isoformat().split('T')[0]
            where_from_date = date_col + " >= " + from_date2.join(['\'', '\''])
        else:
            where_from_date = ''
    else:
        where_from_date = ''

    if isinstance(to_date, str):
        to_date1 = to_datetime(to_date, errors='coerce')
        if isinstance(to_date1, Timestamp):
            to_date2 = to_date1.isoformat().split('T')[0]
            where_to_date = date_col + " <= " + to_date2.join(['\'', '\''])
        else:
            where_to_date = ''
    else:
        where_to_date = ''

    where_stmt.extend([where_from_date, where_to_date])
    where_lst = [i for i in where_stmt if len(i) > 0]
    if len(where_lst) == 0:
        where_lst = None
    return (where_lst)


def sql_ts_agg_stmt(table, groupby_cols, date_col, values_cols, resample_code, period=1, fun='mean', val_round=3,
                    where_lst=None):
    """
    Function to create an SQL statement to pass to an SQL driver to resample a time series table.

    Parameters
    ----------
    table : str
        The SQL table name.
    groupby_cols : str or list of str
        The columns in the SQL table to grouped and returned with the time series data.
    date_col : str
        The date column in the table.
    values_cols : str or list of str
        The column(s) of the value(s) that should be resampled.
    resample_code : str
        The Pandas time series resampling code. e.g. 'D' for day, 'W' for week, 'M' for month, etc.
    period : int
        The number of resampling periods. e.g. period = 2 and resample = 'D' would be to resample the values over a 2 day period.
    fun : str
        The resampling function. i.e. mean, sum, count, min, or max. No median yet...
    val_round : int
        The number of decimals to round the values.
    where_lst : list or None
        A list of where statements to be passed and added to the final SQL statement.

    Returns
    -------
    str
        A full SQL statement that can be passed directly to an SQL connection driver like pymssql through pandas read_sql function.
    """

    if isinstance(groupby_cols, (str, list)):
        groupby_cols_lst = list(groupby_cols)
    else:
        raise TypeError('groupby must be either a str or list of str.')
    if isinstance(values_cols, str):
        values_cols = [values_cols]
    elif not isinstance(values_cols, list):
        raise TypeError('values must be either a str or list of str.')

    pandas_dict = {'D': 'day', 'W': 'week', 'H': 'hour', 'M': 'month', 'Q': 'quarter', 'T': 'minute', 'A': 'year'}
    fun_dict = {'mean': 'avg', 'sum': 'sum', 'count': 'count', 'min': 'min', 'max': 'max'}

    groupby_str = ", ".join(groupby_cols_lst)
    values_lst = ["round(" + fun_dict[fun] + "(" + i + "), " + str(val_round) + ") AS " + i for i in values_cols]
    values_str = ", ".join(values_lst)

    if isinstance(where_lst, list):
        where_stmt = " where " + " and ".join(where_lst)
    else:
        where_stmt = ""

    if isinstance(resample_code, str):
        stmt1 = "SELECT " + groupby_str + ", DATEADD(" + pandas_dict[resample_code] + ", DATEDIFF(" + pandas_dict[
            resample_code] + ", 0, " + date_col + ")/ " + str(period) + " * " + str(
            period) + ", 0) AS " + date_col + ", " + values_str + " FROM " + table + where_stmt + " GROUP BY " + groupby_str + ", DATEADD(" + \
                pandas_dict[resample_code] + ", DATEDIFF(" + pandas_dict[
                    resample_code] + ", 0, " + date_col + ")/ " + str(period) + " * " + str(period) + ", 0)"
    else:
        groupby_cols_lst.extend([date_col])
        groupby_cols_lst.extend(values_cols)
        stmt1 = "SELECT " + ", ".join(groupby_cols_lst) + " FROM " + table + where_stmt

    return (stmt1)


def site_stat_stmt(table, site_col, values_col, fun):
    """
    Function to produce an SQL statement to make a basic summary grouped by a sites column.

    Parameters
    ----------

    table : str
        The database table.
    site_col : str
        The column containing the sites.
    values_col : str
        The column containing the values to be summarised.
    fun : str
        The function to apply.

    Returns
    -------
    str
        SQL statement.
    """
    fun_dict = {'mean': 'avg', 'sum': 'sum', 'count': 'count', 'min': 'min', 'max': 'max'}

    cols_str = ', '.join([site_col, fun_dict[fun] + '(' + values_col + ') as ' + values_col])
    stmt1 = "SELECT " + cols_str + " FROM " + table + " GROUP BY " + site_col
    return (stmt1)


def sql_del_rows_stmt(table, **kwargs):
    """
    Function to create an sql statement to row rows based on where statements.
    """

    where_list = sql_where_stmts(**kwargs)
    stmt1 = "DELETE FROM " + table + " WHERE " + " and ".join(where_list)
    return (stmt1)


def rd_sql_geo(server, database, table, where_lst=None):
    """
    Function to extract the geometry and coordinate system from an SQL geometry field. Returns a shapely geometry object and a proj4 str.

    Parameters
    ----------
    server : str
        The server name. e.g.: 'SQL2012PROD03'
    database : str
        The specific database within the server. e.g.: 'LowFlows'
    table : str
        The specific table within the database. e.g.: 'LowFlowSiteRestrictionDaily'
    where_lst : list
        A list of where statements to be passed and added to the final SQL statement.

    Returns
    -------
    list of shapely geometry objects
        The main output is a list of shapely geometry objects for all queried rows of the SQL table.
    str
        The second output is a proj4 str of the projection system.
    """

    ## Create connection to database
    conn = connect(server, database=database)

    geo_col_stmt = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME=" + "\'" + table + "\'" + " AND DATA_TYPE='geometry'"
    geo_col = str(read_sql(geo_col_stmt, conn).iloc[0, 0])
    geo_srid_stmt = "select distinct " + geo_col + ".STSrid from " + table
    geo_srid = int(read_sql(geo_srid_stmt, conn).iloc[0, 0])
    if where_lst is not None:
        if len(where_lst) > 0:
            stmt2 = "SELECT " + geo_col + ".STGeometryN(1).ToString()" + " FROM " + table + " where " + " and ".join(
                where_lst)
        else:
            stmt2 = "SELECT " + geo_col + ".STGeometryN(1).ToString()" + " FROM " + table
    else:
        stmt2 = "SELECT " + geo_col + ".STGeometryN(1).ToString()" + " FROM " + table
    df2 = read_sql(stmt2, conn)
    df2.columns = ['geometry']
    geo = [loads(x) for x in df2.geometry]
    proj4 = from_epsg_code(geo_srid).to_proj4()

    return (geo, proj4)


#def mssql_logging(log, process):
#    """
#    Function to log to ECan's mssql database logger.
#
#    Parameters
#    ----------
#    log: str
#        Log string.
#    process: str
#        The application process.
#
#    Returns
#    -------
#    None
#    """
#    server = 'SQL2014DEV01'
#    database = 'Ecan_logger'
#    table = 'Log'
#
#    insert_stmt1 = "insert into " + table + "(Date, Application, Thread, Level, Logger, Message) values "
#    dt1 = datetime.now().strftime('%Y-%m-%d %H:%M')
#
#    conn = connect(server, database=database)
#    cursor = conn.cursor()
#
#    log1 = log.replace('\n', ' ').replace('\'', '')
#    rows = "(" + ", ".join([dt1, process, '1', "INFO", "Python", log1]) + ")"
#    insert_stmt2 = insert_stmt1 + rows
#    cursor.execute(insert_stmt2)
#
#    conn.commit()
#
#    #### Close everything!
#    cursor.close()
#    conn.close()















#########################################################
## Archive


