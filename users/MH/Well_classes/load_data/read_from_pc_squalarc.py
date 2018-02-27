"""
Author: matth
Date Created: 7/02/2017 12:35 PM
"""

from __future__ import division

import datetime as dt
import warnings as warn

import pandas as pd

from users.MH.Well_classes.load_data.qual_parameters import qual_parameters, expected_squalarc_unit, convert_qual_units, squalarc_to_standard_names

_KEYS_TO_IGNORE = [
 'datetime',
    'year',
    'month',
'Nitrate-N + Nitrite-N',
'Nitrite-N',
'Total Dissolved Nitrogen',
'Total Nitrogen',
'Total Dissolved Phosphorus',
'Dissolved Oxygen Saturation',
'Conductivity Field',
'Sulphate Sulphur',
'Sulphate Sulphur',
 'Site_ID',
 'Source',
 'Description',
 'SampleID',
 'Status',
 'Date',
 'Time',
 'Project ID',
 'Site_East',
 'Site_North',
 'Sample Comment',
 'Field Notes',
 'Unknown Parameter',
 'Arsenic',
 'Deuterium',
 'Ion Balance',
 'pH',
 'Sum of anions',
 'Sum of cations',
 'Dissolved Organic Carbon',
 'Total Kjeldahl Nitrogen',
 'ORP (mV)',
 'pH Field',
 'Depth to Water',
 'Purging Rate',
 'Purging Time',
 'Sample Appearance',
 'Sample Colour',
 'Sample Odour',
 'Time (NZST)',
 'Field Meter Number',
 'Cadmium',
 'Lead',
 'Copper',
 'Aluminium',
 'Nickel',
 'Boron',
 'Fluoride',
 'Total Suspended Solids',
 'Chromium',
 'Zinc']

# this eventually needs to be integrated with mike's rd_squalarc
def read_from_pc_squalarc_csv_to_dataframe(path, param_list):
    # read in raw data
    # read units and number of lines
    with open(path) as f:
        for i, line in enumerate(f.readlines()):
            if i == 0:
                headers = line.split(',')
            if line.split(',')[0] == '' and ("mg/l" in line.lower() or 'g/m3' in line.lower()):
                units = line.split(',')
                unit_pos = i

    headers = [e.strip() for e in headers]
    units = [e.strip() for e in units]

    unit_dict = {}
    for header, unit in zip(headers, units):
        unit_dict[header] = unit

    data_raw = pd.read_csv(path, nrows=unit_pos, names=headers, skiprows=1, skip_blank_lines=True)
    data_raw = data_raw.loc[pd.notnull(data_raw['Site_ID'])]
    data_raw.index = range(len(data_raw['Site_ID']))
    data_raw['Time'][pd.isnull(data_raw['Time'])] = 1200  # set all NaN time values to noon
    dates = data_raw['Date']
    times = data_raw['Time']

    time_temp = []
    for i, date, time_ in zip(range(0, len(dates)), dates, times):
        try:
            if float(time_) <= 25:  # to handle some weird formats e.g. 10.3 I assume it is a typo
                time_ = float(time_)*100
            time_temp.append(dt.datetime.strptime('{} {:04d}'.format(date, int(float(time_))), "%d-%b-%y %H%M"))
        except:
            try:
                if 'p' in time_.lower():
                    time_ = time_.split(' ')[0] + ' PM'
                elif 'a' in time_.lower():
                    time_ = time_.split(' ')[0] + ' AM'

                time_temp.append(dt.datetime.strptime('{} {}'.format(date, time_), "%d-%b-%y %I:%M:%S %p"))
            except:
                time_temp.append(dt.datetime(1900,1,1)) # a null value, which will be filtered out
                warn.warn('entry {}: {} {} could not be converted to datetime and will be excluded from analysis'.format(i, date, time_))
    years = [e.year for e in time_temp]
    months = [e.month for e in time_temp]
    data_raw['datetime'] = time_temp[:]
    data_raw['year'] = years[:]
    data_raw['month'] = months[:]

    # deal with < and > symbols
    for key in data_raw.keys():
        if key in ['Site_ID', 'Source', 'Description', 'SampleID', 'Status', 'Date', 'Time', 'Project ID', 'Site_East',
                   'Site_North','Sample Comment','Field Notes', 'datetime']:
            continue

        data_raw[key] = data_raw[key].astype(str)
        data_raw[key] = data_raw[key].str.replace('<','-')
        data_raw[key] = data_raw[key].str.replace('>','')
        data_raw[key] = data_raw[key].str.replace(' ','')
        data_raw[key] = pd.to_numeric(data_raw[key],errors='coerce')

    #convert units
    for key in data_raw.keys():
        if key in _KEYS_TO_IGNORE:
            continue
        data_raw[key] = convert_qual_units(data_raw[key],unit_dict[key],expected_squalarc_unit[key])

    # combine values if multiple instances from squalarc dataframe and remove superfluous
    for param in param_list:
        squalarc_names = qual_parameters[param]['squalarc_name']
        if len(squalarc_names) >= 2:
            for i in range(1, len(squalarc_names)):
                idx = pd.isnull(data_raw[squalarc_names[i - 1]])
                data_raw.loc[idx, squalarc_names[0]]= data_raw[squalarc_names[i]][idx]
                data_raw = data_raw.drop(squalarc_names[i], 1)

    #rename series names
    data_raw = data_raw.rename(columns=squalarc_to_standard_names)
    data_raw['proj_id'] = data_raw['proj_id'].str.strip()


    return data_raw


def read_from_pc_squalarc_dataframe(well_num, data_frame, start_date=None, end_date=None):
    """
    reads values from a csv exported from squalarc samples must be in date order rather than grouped by site,
    :param well_num:
    :param data_frame:
    :param start_date: date time object bound selection inclusive
    :param end_date: date time object
    :return:
    """

    #subset by time
    if start_date is not None and end_date is not None:
        data_frame = data_frame.loc[(data_frame['qual_time']>= start_date) & (data_frame['qual_time']<= end_date)]

    #subset by well_num
    data_frame = data_frame.loc[data_frame['Site_ID'] == well_num]

    return data_frame

if __name__ == '__main__':
    path = r"P:\Groundwater\Annual groundwater quality survey 2016\qual_data2016-2017.csv"
    with open(r"P:\Groundwater\Annual groundwater quality survey 2016\annual_survey_well_list.txt") as f:
        well_list = f.readlines()
    well_list = [e.strip() for e in well_list]
    param_list = qual_parameters.keys()
    data_frame = read_from_pc_squalarc_csv_to_dataframe(path, param_list)
    well_data_frame = read_from_pc_squalarc_dataframe('BX24/0347',data_frame)

    print("done")