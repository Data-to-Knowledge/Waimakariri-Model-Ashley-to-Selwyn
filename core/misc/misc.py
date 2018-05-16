# -*- coding: utf-8 -*-
"""
Misc functions for various small procedures.
"""
from __future__ import print_function
from pandas import DataFrame, read_csv, Series, Index, set_option, reset_option
from numpy import ndarray, array, append, isnan, arange
from geopandas import read_file
from re import search, IGNORECASE, findall
import fnmatch, os
from os.path import splitext
from os import path, listdir
from datetime import datetime
import sys


def df_first_valid(df):
    """Get the time index of the first non-na value"""
    def func(x):
        if x.first_valid_index() is None:
            return None
        else:
            return x.first_valid_index()
    df2 = df.apply(func, axis=0)
    return(df2)


def df_last_valid(df):
    """Get the time index of the last non-na value"""
    def func(x):
        if x.last_valid_index() is None:
            return None
        else:
            return x.last_valid_index()
    df2 = df.apply(func, axis=0)
    return(df2)


def printf(x):
    """Print the full rows of a series or dataframe"""
    set_option('display.max_rows', len(x))
    print(x)
    reset_option('display.max_rows')


def rd_dir(data_dir, ext, file_num_names=False, ignore_case=True):
    """
    Function to read a directory of files and create a list of files associated with a spcific file extension. Can also create a list of file numbers from within the file list (e.g. if each file is a station number.)
    """

    if ignore_case:
        files = array([filename for filename in listdir(data_dir) if search('.' + ext + '$', filename, IGNORECASE)])
    else:
        files = array([filename for filename in listdir(data_dir) if search('.' + ext + '$', filename)])

    if file_num_names:
        site_names = array([int(findall("\d+", fi)[0]) for fi in files])
        return([files, site_names])
    else:
        return(files)


def up_branch(df, index_col=1):
    """
    Function to create a dataframe of all the interconnected values looking upstream from specific locations.
    """

    col1 = df.columns[index_col-1]
    index1 = df[col1]
    df2 = df.drop(col1, axis=1)
    catch_set2 = []
    for i in index1:
        catch1 = df2[index1 == i].dropna(axis=1).values[0]
        catch_set1 = catch1
        check1 = index1.isin(catch1)
        while sum(check1) >= 1:
            if sum(check1) > len(catch1):
                print('Index numbering is wrong!')
            catch2 = df2[check1].values.flatten()
            catch3 = catch2[~isnan(catch2)]
            catch_set1 = append(catch_set1, catch3)
            check1 = index1.isin(catch3)
            catch1 = catch3
        catch_set2.append(catch_set1.tolist())

    output1 = DataFrame(catch_set2, index=index1)
    return(output1)


def select_sites(x):
    """
    Function to check for different object types and create an array of values.
    """

    if isinstance(x, ndarray):
        x1 = x.copy()
    elif isinstance(x, (list, tuple)):
        x1 = array(x).copy()
    elif isinstance(x, (Series, Index)):
        x1 = x.values.copy()
    elif isinstance(x, DataFrame):
        x1 = x1.iloc[:, 0].values.copy()
    elif isinstance(x, str):
        if x.endswith('.shp'):
            x1 = read_file(x).copy()
        else:
            x1 = read_csv(x).iloc[:, 0].values.copy()
    elif x is None:
        x1 = x

    return(x1)


def merge_two_dicts(x, y):
    '''
    Given two dicts, merge them into a new dict as a shallow copy.
    '''
    z = x.copy()
    z.update(y)
    return(z)


def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()


def lst_rem_files(path, pattern, rem=False):
    """
    Function to remove all files matching a specific extension.
    """

    # Return all files in dir, and all its subdirectories, ending in pattern
    def gen_files(path, pattern):
       for dirname, subdirs, files in os.walk(path):
          for f in files:
             if f.endswith(pattern):
                yield os.path.join(dirname, f)


    t1 = gen_files(path, pattern)
    # Remove all files matching pattern in the current dir
    if rem:
        for f in t1:
           os.remove(f)
    return(t1)


def grp_mode(df, grp_cols, val_col):
    """
    Groupby mode for Pandas DataFrames.
    """

    df1 = df.groupby(grp_cols)[val_col].value_counts()
    df1.name = 'count'
    df2 = df1.reset_index(val_col).drop('count', axis=1)
    levels = arange(0, len(grp_cols)).tolist()
    df3 = df2.groupby(level=levels)[val_col].first().reset_index()
    return(df3)


def time_switch(x):
    """
    Convenience codes to convert for time text to pandas time codes.
    """
    return {
        'min': 'Min',
        'mins': 'Min',
        'minute': 'Min',
        'minutes': 'Min',
        'hour': 'H',
        'hours': 'H',
        'day': 'D',
        'days': 'D',
        'week': 'W',
        'weeks': 'W',
        'month': 'M',
        'months': 'M',
        'year': 'A',
        'years': 'A',
        'water year': 'A-JUN',
        'water years': 'A-JUN',
    }.get(x, 'A')


#def unarchive_dir(folder, ext='zip', rem_original=False):
#    """
#    Function to unarchive files in all subfolders into those subfolders.
#
#    folder -- The base directory.\n
#    ext -- The archive file extension to be extracted.\n
#    rem_original -- Should the original archive files be removed after extraction?
#    """
#
#    for root, dirs, files in os.walk(folder):
#        for filename in fnmatch.filter(files, '*.' + ext):
#            print(os.path.join(root, filename))
##            pyunpack.zipfile.ZipFile(os.path.join(root, filename)).extractall(root)
#            patoolib.extract_archive(os.path.join(root, filename), outdir=root, interactive=False)
#            if rem_original:
#                os.remove(os.path.join(root, filename))


def save_df(df, path_str, index=True, header=True):
    """
    Function to save a dataframe based on the path_str extension. The path_str must  either end in csv or h5.

    df -- Pandas DataFrame.\n
    path_str -- File path (str).\n
    index -- Should the row index be saved? Only necessary for csv.
    """

    path1 = splitext(path_str)

    if path1[1] in '.h5':
        df.to_hdf(path_str, 'df', mode='w')
    if path1[1] in '.csv':
        df.to_csv(path_str, index=index, header=header)


def get_subdir(a_dir, full_path=False):
    """
    Simple function to get all subdirectories from a directory.
    """
    if full_path:
        return [path.join(a_dir, name) for name in listdir(a_dir) if path.isdir(path.join(a_dir, name))]
    else:
        return [name for name in listdir(a_dir) if path.isdir(path.join(a_dir, name))]


def pytime_to_datetime(pytime):
    """
    Function to convert a PyTime object to a datetime object.
    """

    dt1 = datetime(year=pytime.year, month=pytime.month, day=pytime.day, hour=pytime.hour, minute=pytime.minute)
    return(dt1)


def logging(log_file_path, text):
    """
    Simple function to save/append to a log file.

    Parameters
    ----------
    log_file_path : str
        Path to the log file.
    text : str
        Text str to be added to the log file.

    Returns
    -------
    None
    """

    now1 = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_file_path, 'a+') as file1:
        file1.write(now1 + ' - ' + text + '\n')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



