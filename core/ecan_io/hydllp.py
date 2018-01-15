"""
Functions to read in Hydstra data. Requires a 32bit python environment.
"""

import ctypes
import os
import contextlib
from numpy import array_split, ceil, array
from pandas import to_numeric, to_datetime, concat, DataFrame, Timestamp, merge, to_timedelta, HDFStore, DateOffset
from datetime import date
from core.misc.misc import select_sites, save_df
from core.ecan_io import rd_sql, write_sql
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity


def rd_hydstra(varto, sites=None, data_source='A', from_date=None, to_date=None, from_mod_date=None, to_mod_date=None, interval='day', qual_codes=[30, 20, 10, 11, 21, 18], concat_data=True, export=None):
    """
    Function to read in data from Hydstra's database using HYDLLP. This function extracts all sites with a specific variable code (varto).

    Parameters
    ----------
    varto : int or float
        The hydstra conversion data variable (140.00 is flow).
    sites: list of str
        List of sites to be returned. None includes all sites.
    data_source : str
        Hydstra datasource code (usually 'A').
    from_date: str
        The starting date for the returned data given other constraints.
    to_date: str
        The ending date for the returned data given other constraints.
    from_mod_date: str
        The starting date when the data has been modified.
    to_mod_date: str
        The ending date when the data has been modified.
    interval : str
        The frequency of the output data (year, month, day, hour, minute, second, period). If data_type is 'point', then interval cannot be 'period' (use anything else, it doesn't matter).
    qual_codes : list of int
        The quality codes for output.
    export_path: str
        Path string where the data should be saved, or None to not save the data.

    Return
    ------
    DataFrame
        In long format with site and time as a MultiIndex and data, qual_code, and hydstra_var_code as columns.
    """
    ### Parameters
    device_data_type = {100: 'mean', 140: 'mean', 143: 'mean', 450: 'mean', 110: 'mean', 130: 'mean', 10: 'tot'}

    today1 = date.today()
    dtype_dict = {'Site': 'varchar', 'HydstraCode': 'smallint', 'Time': 'date', 'Value': 'float', 'QualityCode': 'smallint', 'ModDate': 'date'}

    ### Determine the period lengths for all sites and variables
    sites_var_period = hydstra_sites_var_periods(varto=varto, sites=sites, data_source=data_source)
#    sites_list = sites_var_period.site.unique().tolist()
    varto_list = sites_var_period.varto.unique().astype('int32').tolist()

    ### Restrict period ranges - optional
    if isinstance(from_date, str):
        from_date1 = Timestamp(from_date)
        from_date_df = sites_var_period.from_date.apply(lambda x: x if x > from_date1 else from_date1)
        sites_var_period['from_date'] = from_date_df
    if isinstance(to_date, str):
        to_date1 = Timestamp(to_date)
        to_date_df = sites_var_period.to_date.apply(lambda x: x if x > to_date1 else to_date1)
        sites_var_period['to_date'] = to_date_df

    ### Only pull out data according to the modifcation date ranges - optional
    if isinstance(from_mod_date, str):
        sites_block = sites_var_period[sites_var_period.varfrom == sites_var_period.varto]
        varto_block = sites_block.varto.unique().astype('int32').tolist()

        chg1 = hydstra_data_changes(varto_block, sites_block.site.unique(), from_mod_date=from_mod_date, to_mod_date=to_mod_date).drop('to_date', axis=1)
        if 140 in varto_list:
            sites_flow = sites_var_period[(sites_var_period.varfrom != sites_var_period.varto) & (sites_var_period.varto == 140)]
            chg2 = rating_changes(sites_flow.site.unique().tolist(), from_mod_date=from_mod_date, to_mod_date=to_mod_date)
            chg1 = concat([chg1, chg2])

        chg1.rename(columns={'from_date': 'mod_date'}, inplace=True)
        chg3 = merge(sites_var_period, chg1, on=['site', 'varfrom', 'varto'])
        chg4 = chg3[chg3.to_date > chg3.mod_date].copy()
        chg4['from_date'] = chg4['mod_date']
        sites_var_period = chg4.drop('mod_date', axis=1).copy()

    ### Convert datetime to date as str
    sites_var_period2 = sites_var_period.copy()
    sites_var_period2['from_date'] = sites_var_period2['from_date'].dt.date.astype(str)
    sites_var_period2['to_date'] = sites_var_period2['to_date'].dt.date.astype(str)

    site_str_len = sites_var_period2.site.str.len().max()

    if isinstance(export, str):
            if export.endswith('.h5'):
                store = HDFStore(export, mode='a')

    data = DataFrame()
    for tup in sites_var_period2.itertuples(index=False):
        print('Processing site: ' + str(tup.site))
        varto = tup.varto
        data_type = device_data_type[varto]

        df = rd_hydstra_db([tup.site], data_type=data_type, start=tup.from_date, end=tup.to_date, varfrom=tup.varfrom, varto=varto, interval=interval, qual_codes=qual_codes)
        if df.empty:
            continue
        df['HydstraCode'] = varto
        if varto == 143:
            df.loc[:, 'data'] = df.loc[:, 'data'] * 0.001
            df['HydstraCode'] = 140
        ### Make sure the data types are correct
        df.rename(columns={'data': 'Value', 'qual_code': 'QualityCode'}, inplace=True)
        df.index.rename(['Site', 'Time'], inplace=True)
        df.loc[:, 'QualityCode'] = df['QualityCode'].astype('int32')
        df.loc[:, 'HydstraCode'] = df['HydstraCode'].astype('int32')
        df.loc[:, 'ModDate'] = today1
        if isinstance(export, dict):
            df = df.reset_index()
            from_date1 = str(df.Time.min().date())
            to_date1 = str(df.Time.max())
            del_rows_dict = {'where_col': {'Site': [str(tup.site)], 'HydstraCode': [str(df['HydstraCode'][0])]}, 'from_date':from_date1, 'to_date': to_date1, 'date_col': 'Time'}
            write_sql(df, dtype_dict=dtype_dict, del_rows_dict=del_rows_dict, drop_table=False, create_table=False, **export)
        elif isinstance(export, str):
            if export.endswith('.h5'):
                try:
                    store.append(key='var_' + str(varto), value=df, min_itemsize={'site': site_str_len})
                except Exception as err:
                    store.close()
                    raise err
        if concat_data:
            data = concat([data, df])
    if isinstance(export, str):
        store.close()
    if concat_data:
        return data


def rating_changes(sites=None, from_mod_date=None, to_mod_date=None):
    """
    Function to determine flow rating changes during a specified period.

    Parameters
    ----------
    sites: list of str
        List of sites to be returned. None includes all sites.
    from_mod_date: str
        The starting date when the data has been modified.
    to_mod_date: str
        The ending date when the data has been modified.

    Returns
    -------
    DataFrame
        With site, varfrom, varto, and from_date
    """
    ### Parameters
    server= 'SQL2012PROD03'
    database = 'Hydstra'

    table_per = 'RATEPER'
    table_hed = 'RATEHED'
    fields_per = ['STATION', 'VARFROM', 'VARTO', 'SDATE', 'STIME', 'REFTAB', 'PHASE']
    names_per = ['site', 'varfrom', 'varto', 'sdate', 'stime', 'reftab', 'phase']
    fields_hed = ['STATION', 'VARFROM', 'VARTO', 'TABLE', 'RELDATE']
    names_hed = ['site', 'varfrom', 'varto', 'reftab', 'date']

    ### Read data
    if sites is not None:
        if isinstance(sites, list):
            where_col = {'STATION': sites}
        else:
            where_col = None
    else:
            where_col = None

    rate_hed = rd_sql(server, database, table_hed, fields_hed, where_col, rename_cols=names_hed, from_date=from_mod_date, to_date=to_mod_date, date_col='RELDATE')
    rate_hed['site'] = rate_hed['site'].str.strip()

    where_per = {'STATION': rate_hed['site'].astype(str).unique().tolist()}

    rate_per = rd_sql(server, database, table_per, fields_per, where_per, rename_cols=names_per, where_op='OR')
    rate_per['site'] = rate_per['site'].str.strip()
    time1 = to_timedelta(rate_per['stime'] // 100, unit='H') + to_timedelta(rate_per['stime'] % 100, unit='m')
    rate_per['sdate'] = rate_per['sdate'] + time1
    rate_per = rate_per.sort_values(['site', 'sdate']).reset_index(drop=True).drop('stime', axis=1)

    rate_per1 = merge(rate_per, rate_hed[['site', 'reftab']], on=['site', 'reftab'])
    rate_per2 = rate_per1.groupby('site')['sdate'].min().reset_index()
    rate_per2.columns = ['site', 'from_date']

    rate_per2['varfrom'] = 100
    rate_per2['varto'] = 140

    return rate_per2[['site', 'varfrom', 'varto', 'from_date']]


def hydstra_data_changes(varto, sites, data_source='A', from_mod_date=None, to_mod_date=None):
    """
    Function to determine the time series data indexed by sites and variables that have changed between the from_mod_date and to_mod_date. For non-flow rating sites/variables!!!

    Parameters
    ----------
    varto : list of str
        The hydstra conversion data variable (140.00 is flow).
    sites: list of str
        List of sites to be returned.
    data_source : str
        Hydstra datasource code (usually 'A').
    from_mod_date: str
        The starting date when the data has been modified.
    to_mod_date: str
        The ending date when the data has been modified.

    Returns
    -------
    DataFrame
        With site, varfrom, varto, from_date, and to_date
    """
    today1 = Timestamp(date.today())

    ### Get changes for all other parameters
    if isinstance(from_mod_date, str):
        from_mod_date1 = Timestamp(from_mod_date)
        if isinstance(to_mod_date, str):
            to_mod_date1 = Timestamp(to_mod_date)
        else:
            to_mod_date1 = today1
        blocklist = rd_blocklist(sites, [data_source], variables=varto, start_modified=from_mod_date1, end_modified=to_mod_date1)
        if blocklist.empty:
            return blocklist
        else:
            block_grp = blocklist.groupby(['site', 'varto'])
            min_date1 = block_grp['from_mod_date'].min()
            max_date1 = block_grp['to_mod_date'].max()
            min_max_date1 = concat([min_date1, max_date1], axis=1)
            min_max_date1.columns = ['from_date', 'to_date']
            min_max_date2 = min_max_date1.reset_index()
            min_max_date2['varfrom'] = min_max_date2['varto']
            min_max_date3 = min_max_date2[['site', 'varfrom', 'varto', 'from_date', 'to_date']]
            return min_max_date3


def hydstra_sites_var(varto=None, data_source='A', server='SQL2012PROD03', database='Hydstra'):
    """
    Function to extract all of the sites associated with specific varto codes.

    Parameters
    ----------
    varto: list of int or int
        The Hydstra specific variable codes. None equates to all varto's.
    data_source: str
        The Hydstra data source ID. 'A' is archive.

    Returns
    -------
    DataFrame
        With site, data_source, varfrom, and varto
    """
    ### Parameters
    period_tab = 'PERIOD'

    period_cols = ['STATION', 'VARFROM', 'VARIABLE']
    period_names = ['site', 'varfrom', 'varto']

    ## Removals
    rem_dict = {'165131': [140, 140], '69302': [140, 140], '71106': [140, 140], '366425': [140, 140]}

    ### Import

    if varto is None:
        period_where = {'DATASOURCE': data_source}
    elif isinstance(varto, int):
        period_where = {'DATASOURCE': data_source, 'VARIABLE': [varto]}
    elif isinstance(varto, list):
        period_where = {'DATASOURCE': data_source, 'VARIABLE': varto}

    period1 = rd_sql(server, database, period_tab, period_cols, where_col=period_where, rename_cols=period_names)
    period1.loc[:, 'site'] = period1.site.str.strip()

    ### Determine the variables to extract
    period2 = period1[period1.varto.isin(period1.varto.round())].sort_values('site')
    period2 = period2[period2.varto != 101]
    for i in rem_dict:
        period2 = period2[~((period2.site == i) & (period2.varfrom == rem_dict[i][0]) & (period2.varto == rem_dict[i][1]))]

    ### Convert variables to int
    period3 = period2.copy()
    period3['varfrom'] = period3['varfrom'].astype('int32')
    period3['varto'] = period3['varto'].astype('int32')

    ### Return
    return period3


def hydstra_sites_var_periods(varto=None, sites=None, data_source='A', server='SQL2012PROD03', database='Hydstra'):
    """
    Function to determine the record periods for Hydstra sites/variables.

    Parameters
    ----------
    varto : int or float
        The hydstra conversion data variable (140.00 is flow).
    sites: list of str
        List of sites to be returned. None includes all sites.
    data_source : str
        Hydstra datasource code (usually 'A').

    Returns
    -------
    DataFrame
        With site, varfrom, varto, from_date, and to_date.
    """
    sites_var = hydstra_sites_var(varto=varto, data_source=data_source, server=server, database=database)
    if isinstance(sites, list):
        sites_var = sites_var[sites_var.site.isin(sites)]
    sites_list = sites_var.site.unique().tolist()
    varto_list = sites_var.varto.unique().astype('int32').tolist()

    ### Determine the period lengths for all sites and variables
    sites_period = get_site_variable_period(sites_list, data_source)
    sites_period['varfrom'] = sites_period['varto']
    if 140 in varto_list:
        flow_rate_sites = sites_var[(sites_var.varfrom == 100) & (sites_var.varto == 140)]
        wl_sites = sites_period[(sites_period.site.isin(flow_rate_sites.site)) & (sites_period.varto == 100)]
        flow_rate_sites_period = merge(flow_rate_sites, wl_sites[['site', 'from_date', 'to_date']], on='site')

        flow_sites_period = sites_period[sites_period.varto.isin(varto_list)].copy().drop(['var_name', 'units'], axis=1)

        sites_var_period = concat([flow_rate_sites_period, flow_sites_period])
        sites_var_period = sites_var_period[['site', 'varfrom', 'varto', 'from_date', 'to_date']]
    else:
        sites_var_period = sites_period[sites_period.varto.isin(varto_list)].copy().drop(['var_name', 'units'], axis=1)
        sites_var_period = sites_var_period[['site', 'varfrom', 'varto', 'from_date', 'to_date']]

    return sites_var_period.reset_index(drop=True)


def rd_blocklist(sites, datasources=['A'], variables=['100', '10', '110', '140', '130', '143', '450'], start='1900-01-01', end='2100-01-01', start_modified='1900-01-01', end_modified='2100-01-01'):
    """
    Wrapper function to extract info about when data has changed between modification dates.

    Parameters
    ----------
    sites : list, array, one column csv file, or dataframe
        Site numbers.
    datasource : list of str
        Hydstra datasource code (usually ['A']).
    variables : list of int or float
        The hydstra conversion data variable (140.00 is flow).
    start : str
        The start time in the format of '2001-01-01'.
    end : str
        Same formatting as start.
    start_modified: str
        The starting date of the modification.
    end_modified: str
        The ending date of the modification.

    Returns
    -------
    DataFrame
        With site, data_source, varto, from_mod_date, and to_mod_date.
    """
    ### Process sites
    sites1 = select_sites(sites).tolist()

    ### Open connection
    hyd = openHyDb()
    with hyd as h:
        df = h.get_ts_blockinfo(sites1, start=start, end=end, datasources=datasources, variables=variables, start_modified=start_modified, end_modified=end_modified)
    return df


def get_site_variable_period(sites, data_source='A'):
    """

    """
    ### Open connection
    hyd = openHyDb()
    with hyd as h:
        df = h.get_variable_list(sites, data_source)
    return df


def rd_hydstra_db(sites, start=0, end=0, datasource='A', data_type='mean', varfrom=100, varto=140, interval='day', multiplier=1, qual_codes=[30, 20, 10, 11, 21, 18], report_time=None, sites_chunk=20, print_sites=False, export_path=None):
    """
    Wrapper function over hydllp to read in data from Hydstra's database. Must be run in a 32bit python. If either start_time or end_time is not 0, then they both need a date.

    Parameters
    ----------
    sites : list, array, one column csv file, or dataframe
        Site numbers.
    start : str or int of 0
        The start time in the format of either '2001-01-01' or 0 (for all data).
    end : str or int of 0
        Same formatting as start.
    datasource : str
        Hydstra datasource code (usually 'A').
    data_type : str
        mean, maxmin, max, min, start, end, first, last, tot, point, partialtot, or cum.
    varfrom : int or float
        The hydstra source data variable (100.00 is water level).
    varto : int or float
        The hydstra conversion data variable (140.00 is flow).
    interval : str
        The frequency of the output data (year, month, day, hour, minute, second, period). If data_type is 'point', then interval cannot be 'period' (use anything else, it doesn't matter).
    multiplier : int
        interval frequency.
    qual_codes : list of int
        The quality codes in Hydstra for filtering the data.
    sites_chunk : int
        Number of sites to request to hydllp at one time. Do not change unless you understand what it does.

    Return
    ------
    DataFrame
        In long format with site and time as a MultiIndex.
    """

    ### Process sites into workable chunks
    sites1 = select_sites(sites)
    n_chunks = ceil(len(sites1) / float(sites_chunk))
    sites2 = array_split(sites1, n_chunks)

    ### Run instance of hydllp
    data = DataFrame()
    for i in sites2:
        if print_sites:
            print(i)
        ### Open connection
        hyd = openHyDb()
        with hyd as h:
            df = h.get_ts_traces(i, start=start, end=end, datasource=datasource, data_type=data_type, varfrom=varfrom, varto=varto, interval=interval, multiplier=multiplier, qual_codes=qual_codes, report_time=report_time)
        data = concat([data, df])

    if isinstance(export_path, str):
        save_df(data, export_path)

    return data


# Define a context manager generator
# that creates and releases the connection to the hydstra server
@contextlib.contextmanager
def openHyDb(ini_path='Y:/Hydstra/prod/hyd/', dll_path='Y:/Hydstra/prod/hyd/sys/run/', username='', password=''):
    hyd = Hydllp(dll_path=dll_path, ini_path=ini_path)
    try:
        hyd.login(username, password)
        yield hyd
    finally:
        hyd.logout()


# Exception for hydstra related errors
class HydstraError(Exception):
    pass


class HydstraErrorUnknown(HydstraError):
    pass


class Hydllp(object):
    def __init__(self,
                 dll_path,
                 ini_path,
                 hydllp_filename='hydllp.dll',
                 hyaccess_filename='Hyaccess.ini',
                 hyconfig_filename='HYCONFIG.INI'):

        self._dll_path = dll_path
        self._ini_path = ini_path

        self._dll_filename = os.path.join(self._dll_path, hydllp_filename)
        self._hyaccess_filename = os.path.join(self._ini_path, hyaccess_filename)
        self._hyconfig_filename = os.path.join(self._ini_path, hyconfig_filename)

        # See Hydstra Help file
        # According to the HYDLLP doc, the hydll.dll needs to run "in situ" since
        # it needs to reference other files in that directory.
        os.chdir(self._dll_path)

        # According to the HYDLLP doc, the stdcall calling convention is used.
        self._dll = ctypes.WinDLL(self._dll_filename)

        # Hydstra server handle. Unique to each instance.
        self._handle = ctypes.c_int()

        self._logged_in = False

        # ********************************************************************************

    # Start - Define HYDLLP Wrappers
    # ********************************************************************************

    def _decode_error(self, error_code):
        """
        HYDLLP.dll "DecodeError" function.

        Parameters
        ----------
            -error_code, int
                The error code returned by startup_ex
        """

        # Reference the DecodeError dll function
        decode_error_lib = self._dll['DecodeError']
        decode_error_lib.restype = ctypes.c_int

        # string c_type to store the error message
        error_str = ""
        c_error_str = ctypes.c_char_p(error_str)

        # Allocate memory for the return string
        return_str = ctypes.create_string_buffer(" ", 1400)

        # Call "DecodeError"
        err = decode_error_lib(ctypes.c_int(error_code),
                               c_error_str,
                               ctypes.c_int(1023))
        return return_str.value

    def _start_up_ex(self, user, password, hyaccess, hyconfig):
        """
        HYDLLP.dll "StartUpEx" function

        Parameters
        ----------
            -user, str
                Hydstra username

            -password, str
                Hydstra password

            -hyaccess, str
                Fullpath to HYACCESS.INI

            -hyconfig
                Fullpath to HYCONFIG.INI
        """

        startUpEx_lib = self._dll['StartUpEx']
        startUpEx_lib.restype = ctypes.c_int

        # Call the dll function "StartUpEx"
        err = startUpEx_lib(ctypes.c_char_p(user),
                            ctypes.c_char_p(password),
                            ctypes.c_char_p(hyaccess),
                            ctypes.c_char_p(hyconfig),
                            ctypes.byref(self._handle))
        return err

    def _shutdown(self):
        """
        HYDLLP.dll "ShutDown" function

        Parameters
        ----------
            None
        """

        shutdown_lib = self._dll['ShutDown']
        shutdown_lib.restype = ctypes.c_int

        error_code = shutdown_lib(self._handle)

        # Values other than 0 means that an error occured
        if error_code != 0:
            error_msg = self._decode_error(error_code)
            raise HydstraError(error_msg)

    def _json_call(self, request_str, return_str_len):
        """
        HYDLLP.dll "JsonCall" function
        """

        jsonCall_lib = self._dll['JSonCall']
        jsonCall_lib.restype = ctypes.c_int

        # Allocate memory for the return string
        return_str = ctypes.create_string_buffer(" ", return_str_len)

        # c_return_str = ctypes.c_char_p(return_str)

        # Call the dll function "JsonCall"
        err = jsonCall_lib(self._handle,
                           ctypes.c_char_p(request_str),
                           return_str,
                           return_str_len)

        result = return_str.value
        return (result)

        # ********************************************************************************

    # End - Define HYDLLP Wrappers
    # ********************************************************************************

    def login(self, username, password):
        """
        Logs into hydstra using StartUpEx

        Parameters:
        -----------
            -username, str
                Hydstra username

            -passwords, str
                Hydstra password

        """
        error_code = self._start_up_ex(username,
                                       password,
                                       self._hyaccess_filename,
                                       self._hyconfig_filename)

        # Values other than 0 means that an error occured
        if error_code != 0:
            error_msg = self._decode_error(error_code)
            raise HydstraError(error_msg)

        self._logged_in = True

    def logout(self):
        """
        Log out of hydstra

        Parameters:
        ----------
            None
        """
        if self._logged_in:
            self._shutdown()

    def query_by_dict(self, request_dict):
        """
        Sends and receives request to the hydstra server using hydllp.dll.
        """
        import json

        # initial buffer length
        # If it is too small, we can resize, see below
        buffer_len = 1400

        # convert request dict to a json string
        request_json = json.dumps(request_dict)

        # call json_call and convert result to python dictionary
        result_json = self._json_call(request_json, buffer_len)
        result_dict = json.loads(result_json)

        # If the initial buffer is too small, then re-call json_call
        # with the actual buffer length given by the error response
        if result_dict["error_num"] == 200:
            buffer_len = result_dict["buff_required"]
            result_json = self._json_call(request_json, buffer_len)
            result_dict = json.loads(result_json)

        # If error_num is not 0, then an error occured
        if result_dict["error_num"] != 0:
            error_msg = "Error num:{}, {}".format(result_dict['error_num'],
                                                  result_dict['error_msg'])
            raise HydstraError(error_msg)

        # Just in case the result doesn't have a 'return'
        elif 'return' not in result_dict:
            error_msg = "Error code = 0, however no 'return' was found"
            raise HydstraError(error_msg)

        return (result_dict)

    def get_site_list(self, site_list_exp):
        # Generate a request of all the sites
        site_list_req_dict = {"function": "get_site_list",
                              "version": 1,
                              "params": {"site_list": site_list_exp}}

        site_list_result = self.query_by_dict(site_list_req_dict)

        return (site_list_result["return"]["sites"])

    def get_variable_list(self, site_list, data_source):

        # Convert the site list to a comma delimited string of sites
        site_list_str = ",".join([str(site) for site in site_list])

        var_list_request = {"function": "get_variable_list",
                            "version": 1,
                            "params": {"site_list": site_list_str,
                                       "datasource": data_source}}

        var_list_result = self.query_by_dict(var_list_request)
        list1 = var_list_result["return"]["sites"]
        df1 = DataFrame()
        for i in list1:
            site = i['site']
            df_temp = DataFrame(i['variables'])
            df_temp['site'] = site
            df1 = concat([df1, df_temp])

        ## Mangling
        df2 = df1.copy().drop('subdesc', axis=1)
        df2['period_end'] = to_datetime(df2['period_end'], format='%Y%m%d%H%M%S')
        df2['period_start'] = to_datetime(df2['period_start'], format='%Y%m%d%H%M%S')
        df2['site'] = df2['site'].str.strip().astype(str)
        df2['variable'] = df2['variable'].astype(float)
        df2 = df2[df2['variable'].isin(df2['variable'].astype('int32'))]
        df2['variable'] = df2['variable'].astype('int32')
        df3 = df2.drop_duplicates().copy()

        df3.rename(columns={'name': 'var_name', 'period_start': 'from_date', 'period_end': 'to_date', 'variable': 'varto'}, inplace=True)
        df3 = df3[['site', 'varto', 'var_name', 'units', 'from_date', 'to_date']].reset_index(drop=True)

        return df3

    def get_subvar_details(self, site_list, variable):

        # Convert the site list to a comma delimited string of sites
        site_list_str = ",".join([str(site) for site in site_list])

        var_list_request = {"function": "get_subvar_details",
                            "version": 1,
                            "params": {"site_list": site_list_str,
                                       "variable": variable}}

        var_list_result = self.query_by_dict(var_list_request)

        return (var_list_result["return"]["sites"])

    def get_sites_by_datasource(self, data_source):

        # Convert the site list to a comma delimited string of sites
#        data_source_str = ",".join([str(i) for i in data_source])

        var_list_request = {"function": "get_sites_by_datasource",
                            "version": 1,
                            "params": {"datasources": data_source}}

        var_list_result = self.query_by_dict(var_list_request)

        return (var_list_result["return"]["datasources"])

    def get_db_areas(self, area_classes_list):

        db_areas_request = {"function": "get_db_areas",
                            "version": 1,
                            "params": {"area_classes": area_classes_list}}

        db_area_result = self.query_by_dict(db_areas_request)

        return (db_area_result["return"]["sites"])

    def get_ts_blockinfo(self, site_list, datasources=['A'], variables=['100', '10', '110', '140', '130', '143', '450'], start='1900-01-01', end='2100-01-01', start_modified='1900-01-01', end_modified='2100-01-01', fill_gaps=0, auditinfo=0):
        """

        """

        # Convert the site list to a comma delimited string of sites
        sites = select_sites(site_list).astype(str)
        site_list_str = ','.join([str(site) for site in sites])

        ### Datetime conversion
        start = Timestamp(start).strftime('%Y%m%d%H%M%S')
        end = Timestamp(end).strftime('%Y%m%d%H%M%S')
        start_modified = Timestamp(start_modified).strftime('%Y%m%d%H%M%S')
        end_modified = Timestamp(end_modified).strftime('%Y%m%d%H%M%S')

        ### dict request
        ts_blockinfo_request = {"function": "get_ts_blockinfo",
                                "version": 2,
                                "params": {'site_list': site_list_str,
                                           'datasources': datasources,
                                           'variables': variables,
                                           'starttime': start,
                                           'endtime': end,
                                           'start_modified': start_modified,
                                           'end_modified': end_modified
                                           }}

        ts_blockinfo_result = self.query_by_dict(ts_blockinfo_request)
        blocks = ts_blockinfo_result['return']['blocks']
        df1 = DataFrame(blocks)
        if df1.empty:
            return(df1)
        else:
            df1['endtime'] = to_datetime(df1['endtime'], format='%Y%m%d%H%M%S')
            df1['starttime'] = to_datetime(df1['starttime'], format='%Y%m%d%H%M%S')
            df1['variable'] = to_numeric(df1['variable'], errors='coerce', downcast='integer')
            df2 = df1[['site', 'datasource', 'variable', 'starttime', 'endtime']].sort_values(['site', 'variable', 'starttime'])
            df2.rename(columns={'datasource': 'data_source', 'variable': 'varto', 'starttime': 'from_mod_date', 'endtime': 'to_mod_date'}, inplace=True)

            return df2

    def get_ts_traces(self, site_list, start=0, end=0, varfrom=100, varto=140, interval='day', multiplier=1, datasource='A', data_type='mean', qual_codes=[30, 20, 10, 11, 21, 18], report_time=None):
        """

        """

        # Convert the site list to a comma delimited string of sites
        sites = select_sites(site_list).astype(str)
        site_list_str = ','.join([str(site) for site in sites])

        ### Datetime conversion - with dates < 1900
        c1900 = Timestamp('1900-01-01')
        if start != 0:
            start1 = Timestamp(start)
            if start1 > c1900:
                start = start1.strftime('%Y%m%d%H%M%S')
            else:
                start = start1.isoformat(' ').replace('-', '').replace(' ', '').replace(':', '')
        if end != 0:
            end1 = Timestamp(end)
            if end1 > c1900:
                end = end1.strftime('%Y%m%d%H%M%S')
            else:
                end = end1.isoformat(' ').replace('-', '').replace(' ', '').replace(':', '')

        ts_traces_request = {'function': 'get_ts_traces',
                             'version': 2,
                             'params': {'site_list': site_list_str,
                                        'start_time': start,
                                        'end_time': end,
                                        'varfrom': varfrom,
                                        'varto': varto,
                                        'interval': interval,
                                        'datasource': datasource,
                                        'data_type': data_type,
                                        'multiplier': multiplier,
                                        'report_time': report_time}}

        ts_traces_request = self.query_by_dict(ts_traces_request)
        j1 = ts_traces_request['return']['traces']

        ### Convert json to a dataframe
        sites = [str(f['site']) for f in j1]

        out1 = DataFrame()
        for i in range(len(j1)):
            df1 = DataFrame(j1[i]['trace'])
            if not df1.empty:
                df1.rename(columns={'v': 'data', 't': 'time', 'q': 'qual_code'}, inplace=True)
                df1['data'] = to_numeric(df1['data'], errors='coerce')
                df1['time'] = to_datetime(df1['time'], format='%Y%m%d%H%M%S')
                df1['qual_code'] = to_numeric(df1['qual_code'], errors='coerce', downcast='integer')
                df1['site'] = sites[i]
                df2 = df1[df1.qual_code.isin(qual_codes)]
                out1 = concat([out1, df2])

        out2 = out1.set_index(['site', 'time'])[['data', 'qual_code']]

        return out2

# End Class
