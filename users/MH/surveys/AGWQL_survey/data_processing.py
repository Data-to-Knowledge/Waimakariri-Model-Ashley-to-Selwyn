"""
Author: matth
Date Created: 8/02/2017 10:37 AM
"""

from __future__ import division

import warnings as warn

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from users.MH.Well_classes import mav, gv, qual_parameters, read_from_pc_squalarc_csv_to_dataframe, WellQual, get_base_well_data_df
from core.stats.mann_kendall import mann_kendall_test


class AGWQL_survey(object):
    all_param_list = ('samp_id',
                      'proj_id',
                      'tcoli',
                      'ecoli',
                      'NO3_N',
                      'HCO3',
                      'Cl',
                      'SO4',
                      'Ca',
                      'Mg',
                      'Na',
                      'K',
                      'NH4_N',
                      'Fe',
                      'Mn',
                      'SiO2',
                      'drp',
                      'hardness',
                      'pH_field',
                      'pH_lab',
                      'conduct',
                      'DO',
                      'temp')

    trend_param_list = ('samp_id', 'proj_id', 'NO3_N', 'drp')
    trend_year_limit = {'NO3_N':9, 'drp':5}
    all_cwms_zones = (
        'Ashburton',
        'Christchurch - West Melton',
        'Hurunui - Waiau',
        'Kaikoura',
        'Lower Waitaki - South Coastal Canterbury',
        'Orari-Temuka-Opihi-Pareora',
        'Selwyn - Waihora',
        'Upper Waitaki',
        'Waimakariri'
    )
    def __init__(self, year, well_list):
        self._year = year
        self.well_list = np.array(well_list)
        self.well_list.sort()
        self.yearly_wells = {}
        self.trend_wells = {}
        for well in well_list:
            self.yearly_wells[well] = WellQual(well, startdate="{}/08/01".format(year), enddate="{}/02/01".format(year + 1))
            self.trend_wells[well] = WellQual(well, startdate="{}/06/01".format(year - 9), enddate="{}/02/01".format(year + 1))

    def load_data(self, df_yearly=None, df_trend=None):
        # replace with mike's read from squalarc
        print('loading annual quality data')
        all_base_data = get_base_well_data_df(self.well_list)
        for well in self.yearly_wells.values():
            well.load_base_data(df=all_base_data)
            well.load_qual_data(self.all_param_list, data_frame=df_yearly)

        print('loading trend quality data')
        for well in self.trend_wells.values():
            well.load_qual_data(self.trend_param_list, data_frame=df_trend,)
        print('data loaded')

    def clean_data(self):
        print('reducing data to annual GWQ samples')
        self.clean_errors = []

        # for the yearly data
        self._from_bulk_to_annual_sample([self._year], self.yearly_wells, self.all_param_list)

        # for the trend data
        trend_years = range(self._year - 9, self._year + 1)
        self.trend_data = self._from_bulk_to_annual_sample(trend_years,self.trend_wells,self.trend_param_list, data_type='trend')
        print("data has been reduced")

        if len(self.clean_errors) > 0:
            warn.warn('some data was not available, please review self.clean_errors')


    def _from_bulk_to_annual_sample(self, year_list, wells, param_list, data_type = "annual"):
        # initialize dictionary
        missing_well_years = []
        for well in wells.values():
            data = well.qual_data
            out_data = pd.DataFrame(index=year_list, columns=data.keys())

            # create a sample year parameter
            sample_year = []
            for date in data['qual_time']:
                if date.month == 1:
                    sample_year.append(date.year -1)
                else:
                    sample_year.append(date.year)

            data['sample_year'] = np.array(sample_year)

            # subset the annual samples and deal to duplicates
            for year in year_list:
                time_idx = ((data['sample_year']) == year) & ((data['month']<2) | (data['month']>=8))

                pro_id_temp = data['proj_id'][time_idx]

                # try to ID by project ID either GWM_A_N or GWM_A_S
                pro_id_idx = (pro_id_temp.str.contains('GWM')) & (pro_id_temp.str.contains('A'))
                temp = pro_id_temp[pro_id_idx]
                temp_date = data['qual_time'][time_idx][pro_id_idx]

                if len(temp) == 1:
                    for param in param_list:
                        out_data.loc[[year],[param]] = data[param][time_idx][pro_id_idx].iloc[0]
                elif len(temp) == 2 and temp_date.iloc[0].date() == temp_date.iloc[1].date():#dates are the same reduce duplicates
                    for param in param_list:
                        if param in ['samp_id', 'proj_id']:
                            out_data.loc[[year], [param]] = data[param][time_idx][pro_id_idx].iloc[0] #just take the first instance
                            continue
                        temp_param_data = data[param][time_idx][pro_id_idx]
                        if temp_param_data.iloc[0] * temp_param_data.iloc[1] < 0:
                            out_data.loc[[year],[param]] = temp_param_data[(temp_param_data > 0)].iloc[0]
                        else:
                            out_data.loc[[year], [param]] = temp_param_data.mean() # this assumes that LOD values are negative
                elif len(temp) == 2 and data_type == 'trend':
                    for param in param_list:
                        if param in ['samp_id', 'proj_id']:
                            out_data.loc[[year], [param]] = data[param][time_idx][pro_id_idx].iloc[
                                0]  # just take the first instance
                            continue
                        temp_param_data = data[param][time_idx][pro_id_idx]
                        temp_param_data_pos = temp_param_data[temp_param_data > 0]
                        if len(temp_param_data_pos) == 0:
                            out_data.loc[[year], [param]] = temp_param_data.iloc[0]
                        else:
                            out_data.loc[[year], [param]] = temp_param_data_pos.mean()
                elif len(temp) == 2:
                    missing_well_years.append(("2 entries 2 date",well.well_num,year))
                elif len(temp) > 2:
                    missing_well_years.append(('>2 entries',well.well_num,year))
                elif len(temp) == 0:
                    if len(pro_id_temp) == 0:
                        continue # if no data for the year do not record as a "clean error"
                    if data_type == 'trend':
                        time_idx2 = ((data['sample_year']) == year) & ((data['month'] < 2) | (data['month'] >= 9))
                        pro_id_temp2 = data['proj_id'][time_idx2]
                        # pull data from all GWM samples in the date range where all are <lod use lod otherwise mean
                        # of detection values.  This assumes LOD values are negative
                        pro_id_idx = (pro_id_temp2.str.contains('GWM',na=False)) | (pro_id_temp2.str.contains('GWI',na=False)) # check this
                        temp = pro_id_temp2[pro_id_idx]
                        if len(temp) == 0:
                            missing_well_years.append(('no entries trend', well.well_num, year))
                        else:
                            for param in param_list:
                                if param in ['samp_id', 'proj_id']:
                                    out_data.loc[[year], [param]] = data[param][time_idx2][pro_id_idx].iloc[
                                        0]  # just take the first instance
                                    continue
                                temp_param_data = data[param][time_idx2][pro_id_idx]
                                temp_param_data_pos = temp_param_data[temp_param_data > 0]
                                if len(temp_param_data_pos) == 0:
                                    out_data.loc[[year], [param]] = temp_param_data.iloc[0]
                                else:
                                    out_data.loc[[year], [param]] = temp_param_data_pos.mean()
                    else:
                        missing_well_years.append(('no entries',well.well_num,year))
                else:
                    missing_well_years.append((3,well.well_num,year))
            well.qual_data = out_data
        self.clean_errors.extend( missing_well_years)

    def calc_trends(self):
        print("calculating mann-kendall trends")
        for well in self.well_list:
            for param in self.trend_param_list:
                if param in ['proj_id', 'samp_id']:
                    continue

                if all(pd.isnull(self.yearly_wells[well].get_param(param))):
                    self.yearly_wells[well].add_param_qual("{}_trend".format(param), np.nan)
                    self.yearly_wells[well].add_param_qual("{}_trend_num_samp".format(param), np.nan)
                    continue

                data = np.array(self.trend_wells[well].get_param(param))
                data.astype(float)
                data = data[pd.notnull(data)]
                data[data<0] *= -0.5
                data_len = (len(data))
                if len(data) < self.trend_year_limit[param]:
                    #can run with as little as four, but 6 or more are recommended
                    trend = 'ne_data'
                else:
                    trend, h, p, z = mann_kendall_test(data, alpha=0.1)
                self.yearly_wells[well].add_param_qual("{}_trend".format(param), trend)
                self.yearly_wells[well].add_param_qual("{}_trend_num_samp".format(param), data_len)

    def save_output_tables(self, out_dir):
        """
        note all_output.csv the LOD values are negative
        :param out_dir:
        :return:
        """
        print('saving output data')
        # spread sheet for GIS integration and overview data
        output_keys = ['lat', 'lon', 'cwms_zone', 'depth']  # Well_base data of interest
        temp_param_list = list(self.trend_param_list)
        temp_param_list.remove('samp_id')
        temp_param_list.remove('proj_id')
        output_keys.extend(["{}_trend".format(e) for e in temp_param_list])
        output_keys.extend(["{}_trend_num_samp".format(e) for e in temp_param_list])
        output_keys.extend(self.all_param_list)
        all_output = pd.DataFrame(index=self.well_list, columns=output_keys,)
        for well in self.well_list:
            for key in output_keys:
                temp = self.yearly_wells[well].get_param(key)
                if isinstance(temp,pd.Series):
                    temp = temp.iloc[0]
                all_output.loc[well,key] = temp
        all_output.to_csv(out_dir + 'all_output.csv')

        # save supporting info
        self._save_supporting_info(out_dir,all_output)

        # regional survey stats (parameter, median and range)
        reg_sum_columns = ['Water Quality Parameters', 'Units', 'Median', 'Range']
        reg_sum = pd.DataFrame(index=range(len(self.all_param_list)), columns=reg_sum_columns)
        for i, param in enumerate(self.all_param_list):
            if param in ['samp_id', 'proj_id']:
                continue
            reg_sum.loc[i,'Water Quality Parameters'] = qual_parameters[param]['standard_name']
            reg_sum.loc[i,'Units'] = qual_parameters[param]['unit']

            values = np.array(all_output[param]).astype(float)


            # median and range
            if np.nanmin(values) < 0:
                lod = np.nanmin(values) * -1
            else:
                lod = 0
            values[values < 0] *= -0.5
            med_rang = {}
            med_rang['Median'] = np.nanmedian(values)
            med_rang['Min'] = np.nanmin(values)
            med_rang['Max'] = np.nanmax(values)
            numdig = qual_parameters[param]['sig_fig']
            for key in med_rang:
                if med_rang[key] < lod:
                    med_rang[key] = round(med_rang[key]*2, numdig)
                    if numdig == 0:
                        med_rang[key] = int(med_rang[key])
                    med_rang[key] = "<{}".format(med_rang[key])
                else:
                    med_rang[key] = round(med_rang[key], numdig)
                    if numdig == 0:
                        med_rang[key] = int(med_rang[key])
                    med_rang[key] = "{}".format(med_rang[key])

            reg_sum.loc[i, 'Median'] = med_rang['Median']
            reg_sum.loc[i, 'Range'] = "{} to {}".format(med_rang['Min'], med_rang['Max'])

        reg_sum.to_csv(out_dir + "regional_summary.csv")

        # comparison to drinking water quality for total and by CWMSZone

        #MAV
        mav_params = ['ecoli', 'NO3_N', 'Mn']
        MAV = pd.DataFrame(index=range(len(mav_params)+1),columns=['parameter', 'all Canterbury'] + list(self.all_cwms_zones))
        for i,param in enumerate(mav_params):
            i = i+1 # offset so that I can have the well count first
            MAV.loc[i,'parameter'] = qual_parameters[param]['standard_name']
            tempdata = all_output[pd.notnull(all_output[param])]
            # all canterbury
            count = (tempdata[param] > mav [param]).sum()
            MAV.loc[i,'all Canterbury'] = count

            # by zone
            for cwms in self.all_cwms_zones:
                tempdata2 = tempdata[tempdata['cwms_zone'] == cwms]
                count = (tempdata2[param] > mav[param]).sum()
                MAV.loc[i,cwms] = count

        # count wells sampled
        MAV.loc[0,'parameter'] = 'well count'
        tempdata = all_output.dropna(how='all', subset=self.all_param_list)
        #all canterbury
        MAV.loc[0,'all Canterbury'] = len(tempdata['samp_id'])
        for cwms in self.all_cwms_zones:
            tempdata2 = tempdata[tempdata['cwms_zone'] == cwms]
            count = len(tempdata2['samp_id'])
            MAV.loc[0, cwms] = count

        #GV
        gv_params = ['NH4_N', 'Cl', 'hardness', 'Fe', 'Mn', 'pH_lab', 'Na', 'SO4']
        GV = pd.DataFrame(index=range(len(gv_params)),columns=['parameter', 'all Canterbury'] + list(self.all_cwms_zones))
        for i,param in enumerate(gv_params):
            GV.loc[i, 'parameter'] = qual_parameters[param]['standard_name']
            tempdata = all_output[pd.notnull(all_output[param])]
            # all canterbury
            if param == 'pH_lab':
                count = ((tempdata[param] < gv[param][0]) | (tempdata[param] > gv[param][1])).sum()
            else:
                count = (tempdata[param] > gv[param]).sum()

            GV.loc[i, 'all Canterbury'] = count

            # by zone
            for cwms in self.all_cwms_zones:
                tempdata2 = tempdata[tempdata['cwms_zone'] == cwms]
                if param == 'pH_lab':
                    count = ((tempdata2[param] < gv[param][0]) | (tempdata2[param] > gv[param][1])).sum()
                else:
                    count = (tempdata2[param] > gv[param]).sum()
                GV.loc[i,cwms] = count

        MAV.to_csv(out_dir + 'MAV_counts.csv')
        GV.to_csv(out_dir + 'GV_counts.csv')

    def _save_supporting_info(self, outdir, all_output):
        # Text file with:
        with open(outdir+'supportinginfo.txt', 'w') as f:
            #  of wells <= 3.0 NO3 (also percentage)
            num_low_n = len(all_output.NO3_N[all_output.NO3_N <=3])
            num_well_samp_n = pd.notnull(all_output.NO3_N).sum()
            f.write('{} wells ({}% of sampled had NO3 less than 3mg/L\n'.format(num_low_n,num_low_n/num_well_samp_n*100))
            # of wells with enough data to do trends
            num_n_trends = ((all_output.NO3_N_trend == 'increasing') | (all_output.NO3_N_trend == 'no trend') |
                            (all_output.NO3_N_trend == 'decreasing')).sum()
            num_p_trends = ((all_output.drp_trend == 'increasing') | (all_output.drp_trend == 'no trend') |
                            (all_output.drp_trend == 'decreasing')).sum()
            f.write('{} n and {}p had enough data to calculate trends\n'.format(num_n_trends,num_p_trends))
            # wells over GV for either mn or fe and both
            fe_gv = ((all_output.Fe > gv['Fe']) & (all_output.Mn <= gv['Mn'])).sum()
            mn_gv = ((all_output.Fe < gv['Fe']) & (all_output.Mn > gv['Mn'])).sum()
            both_gv = ((all_output.Fe > gv['Fe']) & (all_output.Mn > gv['Mn'])).sum()
            f.write('{} wells over GV for Fe, {} wells over GV for Mn, {} wells over GV for both\n'.format(fe_gv,mn_gv,both_gv))
            # ecoli detections <20m depth, number with detections % wells over 50m that had ecoli detections
            e_20 = ((all_output.depth < 20) & (all_output.ecoli > 0)).sum()
            e_50 = ((all_output.depth > 50) & (all_output.ecoli > 0)).sum() / ((all_output.depth > 50) & (pd.notnull(all_output.ecoli))).sum() * 100
            e_over_20 = ((all_output.depth >= 20) & (all_output.ecoli > 0)).sum() / ((all_output.depth >= 20) & (pd.notnull(all_output.ecoli))).sum() * 100
            f.write('ecoli: {} detections <20m, {} percent in wells  >50m, {} percent in wells >=20\n'.format(e_20, e_50, e_over_20))

        # tables
        # trends by zone
        data_columns = ['drp_increasing', 'drp_no trend', 'drp_decreasing', 'NO3_N_increasing', 'NO3_N_no trend', 'NO3_N_decreasing']
        trend_by_zone = pd.DataFrame(index=['all canterbury'] + list(self.all_cwms_zones), columns=data_columns)
        for param in self.trend_param_list:
            if param in ['samp_id', 'proj_id']:
                continue
            for val in ['increasing','no trend', 'decreasing']:
                key = param + '_trend'
                for zone in ['all canterbury'] + list(self.all_cwms_zones):
                    if zone == 'all canterbury':
                        count = (all_output[key] == val).sum()
                    else:
                        count = (all_output[key][all_output.cwms_zone == zone] == val).sum()
                    trend_by_zone.loc[zone,param+ '_'+ val] = count
        trend_by_zone.to_csv(outdir+'trend_summary.csv')

        # diginosise plots
        # trends vs depth
        for param in self.trend_param_list:
            if param in ['samp_id', 'proj_id']:
                continue

            trend = all_output[param + '_trend']
            trend = trend.replace('ne_data', np.nan)
            trend, vals = pd.factorize(trend)
            trend = trend.astype(float)
            trend[trend < 0] = np.nan
            depth = np.array(all_output.depth)
            fig = plt.figure()
            plt.plot(depth, trend, 'ro')
            plt.xlabel('depth')
            plt.ylabel('{}trend'.format(param))
            plt.yticks(range(len(vals)),vals)
            fig.savefig(outdir+param+'trend.png')

        # DRP vs depth
        depth = np.array(all_output.depth)
        drp = np.array(all_output.drp)
        fig = plt.figure()
        plt.plot(depth, drp,'ro')
        plt.xlabel('depth')
        plt.ylabel('drp')
        fig.savefig(outdir+'drpvsdepth.png')


if __name__ == '__main__':
    #input
    #well_list
    well_list_path = r"P:\Groundwater\Annual groundwater quality survey 2016_true_copy\Well_list_2017-03-17.txt"

    #squalarc csvs
    trend_path = r"P:\Groundwater\Annual groundwater quality survey 2016_true_copy\qual_data_trend_2004-2017_3.csv"
    yearly_path = r"P:\Groundwater\Annual groundwater quality survey 2016_true_copy\qual_data2016-2017_unselected_3.csv"


    with open(well_list_path) as f:
        well_list = f.readlines()
    well_list = [e.strip() for e in well_list]
    test = AGWQL_survey(2016,well_list)
    df_trend = read_from_pc_squalarc_csv_to_dataframe(trend_path,test.trend_param_list)
    df_yearly = read_from_pc_squalarc_csv_to_dataframe(yearly_path, test.all_param_list)
    test.load_data(df_trend=df_trend,df_yearly=df_yearly)

    for well in test.trend_wells.values():
        test2 = set([well.well_num]) != set(well.qual_data['Site_ID']) and not set(well.qual_data['Site_ID']) == set()
        if test2:
            warn.warn("wrong well {}: {}".format(well.well_num,set(well.qual_data['Site_ID'])))
    test.clean_data()
    problem_wells = [e[1] for e in test.clean_errors]
    problem_wells_set = list(set(problem_wells))
    problems = []
    for well in problem_wells_set:
        temp = []
        for problem in test.clean_errors:
            if well in problem:
                temp.append((problem[0],problem[2]))
        problems.append((well,temp))

    with open("P:/Groundwater/Annual groundwater quality survey 2016/data/AGWQL.txt", mode='w') as f:
        [f.write(str(e)+'\n') for e in test.clean_errors]

    with open("P:/Groundwater/Annual groundwater quality survey 2016/data/AGWQL_sorted.txt", mode='w') as f:
        [f.write(str(e)+'\n') for e in problems]

    test.calc_trends()
    test.save_output_tables("P:/Groundwater/Annual groundwater quality survey 2016/data/")