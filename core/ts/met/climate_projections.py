# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:02:29 2017

@author: MichaelEK
"""
from xarray import open_dataset, concat
from os import path, walk, makedirs
from pandas import to_datetime


def proc_waimak_nc(base_path, out_path):
    """
    Function to process TopNet output data that come as multiple nc files of different time periods.

    base_path -- The base path where many nc files and subfolders lay.\n
    out_path -- The output directory for the processed data.\n
    """

    #### Walk through the files and folders
    for root, dirs, files in walk(base_path):
        files2 = [i for i in files if i.endswith('.nc')]
        file_paths1 = [path.join(root, i) for i in files2]
        if len(file_paths1) > 0:
            new_base_path = root.replace(base_path, out_path)
            if not path.exists(new_base_path):
                makedirs(new_base_path)
            ### Load in each file and process each
            mtypes = ['flow', 'drainge', 'ovstn_q', 'soilevp', 'soilh2o', 'aprecip']
            for i in file_paths1:
                tds1 = open_dataset(i)

                ## Reorganise variables and coordinates
                tds1['nrch'] = tds1['rchid']
                if 'river_flow_rate_mod' in tds1.data_vars.keys():
                    tds1.rename({'river_flow_rate_mod': 'flow'}, inplace=True)
                if 'mod_flow' in tds1.data_vars.keys():
                    tds1.rename({'mod_flow': 'flow'}, inplace=True)
                tds2 = tds1[mtypes]
                tds2a = tds2.drop(['basin_lat', 'basin_lon', 'end_lat', 'end_lon'])
                tds3 = tds2a.squeeze('nens')
                time1 = to_datetime(to_datetime(tds3['time'].data).date)
                tds3['time'] = time1

                ## Data integrety checks
                freq1 = time1.inferred_freq
                if freq1 != 'D':
                    raise ValueError('Time dimension is not regular!')
                if 'topds1' in locals():
                    force1 = topds1.attrs['forcing']
                    force2 = tds3.attrs['forcing']
                    if force1 != force2:
                        print(path.split(i)[1] + ' has a different forcing compared to previous files for this model')

#                print(path.split(i)[1], tds1.data_vars)

                ## Combine data
                if 'topds1' in locals():
                    topds1 = concat([topds1, tds3], 'time')
                else:
                    topds1 = tds3

            topds2 = topds1.copy()
            del(topds1)

            ### Extract necessary attributes
            model_name = path.split(new_base_path)[1]
            rcp_name = path.split(path.split(new_base_path)[0])[1]
            dates = to_datetime(topds2.time.data)
            start_date = str(dates[0].year)
            end_date = str(dates[-1].year)
#            base_file_name = '_'.join([start_str, rcp_name, model_name]) + '.nc'

            ### Save the new data
            for j in mtypes:
                file_name = '_'.join([j, model_name, rcp_name, start_date, end_date]) + '.nc'
                topds2[[j]].to_netcdf(path.join(new_base_path, file_name))
                print(file_name)

            tds1.close()
            topds2.close()


def proc_cant_nc(base_path, out_path):
    """
    Function to process TopNet output data that come as multiple nc files of different time periods.

    base_path -- The base path where many nc files and subfolders lay.\n
    out_path -- The output directory for the processed data.\n
    """

    mtypes = ['mod_flow', 'drainge', 'canevap', 'instn_q', 'ovstn_q', 'soilevp', 'soilh2o', 'aprecip', 'potevap', 'soilevp', 'soilh2o', 'zbarh2o']

    #### Walk through the files and folders
    for root, dirs, files in walk(base_path):
        files2 = [l for l in files if l.endswith('.nc')]
        file_paths1 = [path.join(root, p) for p in files2]
        if len(file_paths1) > 0:

            ### Make new paths
            new_base_path = root.replace(base_path, out_path)
            if not path.exists(new_base_path):
                makedirs(new_base_path)

            ### Extract necessary attributes
            model_name = path.split(new_base_path)[1]
            rcp_name = path.split(path.split(new_base_path)[0])[1]

            ### Check to see how the files are organised
            for i in mtypes:
                file1 = [j for j in file_paths1 if i in j]
                n_file = len(file1)
                if n_file == 0:
                    raise ValueError('No files for ' + i)
                elif n_file > 1:
                    for f in file1:
                        ttds1 = open_dataset(f)
                        ttds1['nrch'] = ttds1['rchid']
                        ttds2 = ttds1[[i]]
                        ttds2 = ttds2.drop(['end_lat', 'end_lon'])
                        ttds2 = ttds2.squeeze('nens')
                        if 'tds3' in locals():
                            tds3 = concat([tds3, ttds2], 'time')
                        else:
                            tds3 = ttds2

                elif n_file == 1:
                    tds1 = open_dataset(file1[0])
                    tds1['nrch'] = tds1['rchid']
                    tds1 = tds1.drop(['end_lat', 'end_lon'])
                    tds2 = tds1[[i]]
                    tds3 = tds2.squeeze('nens')

                ## Process nc file
                if i == 'mod_flow':
                    tds3.rename({'mod_flow': 'flow'}, inplace=True)
                    i = 'flow'
                time1 = to_datetime(to_datetime(tds3['time'].data).date)
                tds3['time'] = time1
                tds3.attrs['mtype'] = i
                tds3.attrs['rcp'] = rcp_name
                tds3.attrs['model'] = model_name

                ## Data integrety checks
                freq1 = time1.inferred_freq
                if freq1 != 'D':
                    raise ValueError('Time dimension is not regular!')

                ## Extract other parameters for renaming
                dates = to_datetime(tds1.time.data)
                start_date = str(dates[0].year)
                end_date = str(dates[-1].year)

                ## Copy file to new location
                file_name = '_'.join([i, model_name, rcp_name, start_date, end_date, 'topnet_canterbury']) + '.nc'
                new_file = path.join(new_base_path, file_name)

                ## DETERMINE AND ADD IN PROPER COMPRESSION SETTINGS!!!

                tds3.to_netcdf(new_file)
                tds1.close()
                tds3.close()
                del(tds3)



            topds2 = topds1.copy()
            del(topds1)

            ### Extract necessary attributes
            model_name = path.split(new_base_path)[1]
            rcp_name = path.split(path.split(new_base_path)[0])[1]
            dates = to_datetime(topds2.time.data)
            start_date = str(dates[0].year)
            end_date = str(dates[-1].year)
#            base_file_name = '_'.join([start_str, rcp_name, model_name]) + '.nc'

            ### Save the new data
            for j in mtypes:
                file_name = '_'.join([j, model_name, rcp_name, start_date, end_date]) + '.nc'
                topds2[[j]].to_netcdf(path.join(new_base_path, file_name))
                print(file_name)

            tds1.close()
            topds2.close()


















