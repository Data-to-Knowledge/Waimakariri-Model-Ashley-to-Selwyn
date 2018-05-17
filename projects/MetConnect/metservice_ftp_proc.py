# -*- coding: utf-8 -*-
"""
Created on Fri Jul 07 09:09:11 2017

@author: MichaelEK

Script to process the metservice netcdf files.
"""
import sys
from os import path, getcwd
import sqlalchemy

try:

    ######################################
    ### Change pythonpath and add init

    base_py_path = r'D:\Executables\PythonScripts'
    proj_name = 'MetConnect'
    py_path = path.join(base_py_path, proj_name)

    sys.path.append(py_path)

    #init1 = path.join(py_path, '__init__.py')
    #fh = open(init1,'a+')

    #######################################
    ### Start the work

    from configparser import ConfigParser
    from pandas import to_datetime
    from pandas.io.sql import DatabaseError
    from core.ts.met.metservice import proc_metservice_nc, MetS_nc_to_df, metconnect_id_loc
    from core.spatial.vector import sel_sites_poly
    from core.spatial.raster import point_interp_ts
    from core.misc import rd_dir, logging
    from core.ecan_io import write_sql, rd_sql


    def create_engine(db_type, server, database, username=None, password=None):
        """
        Function to create an sqlalchemy engine.

        Parameters
        ----------
        db_type : str
            The type of database to connect to. Options include mssql, postgresql, oracle, mysql, and sqlite.
        server : str
            The server name. e.g.: 'SQL2012PROD03'
        database : str
            The specific database within the server. e.g.: 'LowFlows'
        username : str or None
            Either the username or None when not needed.
        password : str or None
            Either the password or None when not needed.

        Returns
        -------
        sqlalchemy engine

        Notes
        -----
        If pymssql is installed, create_eng will use the package instead of pyodbc.
        """
        if isinstance(username, str):
            up = username
            if isinstance(password, str):
                up = up + ':' + password
            up = up + '@'
        else:
            up = ''
        if db_type == 'mssql':
            try:
                import pymssql
                eng_str = 'mssql+pymssql://' + up + server + '/' + database
                engine = sqlalchemy.create_engine(eng_str)
            except:
                driver1 = '?driver=ODBC+Driver+13+for+SQL+Server'
                eng_str = 'mssql+pyodbc://' + up + server + '/' + database + driver1
                engine = sqlalchemy.create_engine(eng_str)
                try:
                    engine.connect()
                except:
                    driver2 = '?driver=ODBC+Driver+11+for+SQL+Server'
                    eng_str = 'mssql+pyodbc://' + up + server + '/' + database + driver2
                    engine = sqlalchemy.create_engine(eng_str)
                    try:
                        engine.connect()
                    except:
                        driver2 = '?driver=SQL+Server+Native+Client+11.0'
                        eng_str = 'mssql+pyodbc://' + up + server + '/' + database + driver2
                        engine = sqlalchemy.create_engine(eng_str)
                        try:
                            engine.connect()
                        except:
                            print('Install a proper ODBC mssql driver')
        elif db_type == 'postgresql':
            eng_str = 'postgresql://' + up + server + '/' + database
            engine = sqlalchemy.create_engine(eng_str)
        elif db_type == 'oracle':
            eng_str = 'oracle://' + up + server + '/' + database
            engine = sqlalchemy.create_engine(eng_str)
        elif db_type == 'mysql':
            eng_str = 'mysql+mysqldb://' + up + server + '/' + database
            engine = sqlalchemy.create_engine(eng_str)
        elif db_type == 'sqlite':
            engine = sqlalchemy.create_engine('sqlite:///:memory:')

        return engine

    ##########################################
    ### Load in ini parameters

    py_dir = path.realpath(path.join(getcwd(), path.dirname(__file__)))

    ini1 = ConfigParser()
    ini1.read([path.join(py_dir, path.splitext(__file__)[0] + '.ini')])

    nc_dir = str(ini1.get('Input', 'nc_dir'))
    server = str(ini1.get('Input', 'server'))
    data_table = str(ini1.get('Output', 'data_table'))
    nc_output_dir = str(ini1.get('Output', 'nc_output_dir'))
    log_path = str(ini1.get('Output', 'log_path'))

    ## Processing parameters

    database = 'MetConnect'
    mc_site_table = 'RainFallPredictionSitesGrid'
    mc_cols = ['MetConnectID', 'SiteString', 'TidedaID']

    point_site_col = 'MetConnectID'
    time_col = 'time'
    y_col = 'y'
    x_col = 'x'
    data_col = 'precip'
    buffer_dis = 16000
    digits = 2

    ## Output parameters
    create_table = False

    rename_dict = {'site': 'MetConnectID', 'model_date': 'PredictionDateTime', 'time': 'ReadingDateTime', 'precip': 'HourlyRainfall'}

    dtype_dict = {'MetConnectID': 'INT', 'PredictionDateTime': 'DATETIME', 'ReadingDateTime': 'DATETIME', 'HourlyRainfall': 'decimal(18,2)'}

    #######################################
    ### Read in site locations
    print('Read in site locations')

    points = metconnect_id_loc(mc_server=server, mc_db=database, mc_site_table=mc_site_table, mc_cols=mc_cols, gis_server='SQL2012PROD04')

    ########################################
    ### Check the model of the day! And see if we already have a record.
    print('Check the model of the day! And see if we already have a record.')

    txt = rd_dir(nc_dir, 'txt')[-1]
    model1 = open(path.join(nc_dir, txt), 'r').readline()

    nc1 = rd_dir(nc_dir, 'nc')

    nc_model = [i for i in nc1 if model1 in i][-1]
    model_date = to_datetime(nc_model[-16:-6], format='%Y%m%d%H')

    last_predict_date = rd_sql(server, database, stmt = str("select max(PredictionDateTime) from " + data_table)).loc[0][0]

    if model_date == last_predict_date:
        logging(log_path, str(last_predict_date) + ' was the last model date in the database table. It is the lastest date on the ftp site. Nothing was updated.')
        sys.exit('Already have the data for this forecast...closing...')
    else:
        print("We don't have the newest data in the database table...processing data now...")


    #########################################
    ### Process data

    ## preprocess the nc file to save it as a proper nc file
    new_nc = proc_metservice_nc(path.join(nc_dir, nc_model), export_dir=nc_output_dir)

    ## Extract the data from the new nc file
    precip, sites, start_date = MetS_nc_to_df(new_nc)

    ## Select the precip data within the buffer area of the points shp
    #points = read_file(point_shp)
    points_poly1 = points.to_crs(sites.crs).buffer(buffer_dis)
    sites2 = sel_sites_poly(sites, points_poly1)

    precip2 = precip[precip.site.isin(sites2.site)]

    ## Interpolate the data at the shp points
    data = point_interp_ts(precip2, time_col, x_col, y_col, data_col, points, point_site_col, sites.crs, to_crs=None, interp_fun='cubic', agg_ts_fun=None, period=None, digits=2)

    ########################################
    ### Output data
    print('Saving data')

    ## Reformat
    data.loc[:, 'model_date'] = start_date
    data1 = data[['site', 'model_date', 'time', 'precip']].copy()
    data2 = data1.rename(columns=rename_dict)

    ##Output
    ### Prepare the engine
    engine = create_engine('mssql', server, database)

    ### Save to mssql table
    data2.to_sql(name=data_table, con=engine, if_exists='append', chunksize=1000, index=False)

    print('Success!')
    logging(log_path, str(last_predict_date) + ' was the last model date in the database table. ' + str(model_date) + ' is the new model date. ' + str(len(points)) + ' sites were updated.')

except Exception as err:
    err1 = err
    print(err1)
    logging(log_path, 'Run failed with error: ' + str(err1))
    raise err1
