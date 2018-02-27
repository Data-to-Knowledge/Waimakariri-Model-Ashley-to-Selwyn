# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 14:09:36 2017

@author: KateSt
"""

from core.classes.hydro import hydro, all_mtypes
import calendar
from pandas import DateOffset
from numpy import logical_and


site = 64304 #site number
qual_codes = [10, 18, 20, 30, 50]

fr_flow = 0.7 #full restriction flow
pr_flow = 0.910 #partial restirction flow

sites1 = [site] #change site to list so it can be passed to hydro class function
h1 = hydro().get_data(mtypes='flow', sites=sites1, qual_codes=qual_codes) #use hydroclass function to extract flow data

df1 = h1.sel_ts(mtypes='flow', pivot=True) #reformat hydroclass data to dataframe
df2 = df1

findgaps = df2.dropna() #drop days with no flow from the dataset
findgaps['start_gap'] = findgaps.index #find the final day with flow before a gap
findgaps['end_gap'] = findgaps['start_gap'].shift(periods=-1) + DateOffset(-1)
 #find the day preceeding the next day with flow data (this will be the same as start_gap if no gap exists)
findgaps['diff'] = abs(findgaps['start_gap'] - findgaps['end_gap'])
#Find how long the gaps are (will be '0 days' if there is no gap)
gaps = findgaps.loc[findgaps['diff'] > '0 days'] #find actual gaps

gap_number = len(gaps) #find the number of gaps in the record

for i in range(gap_number): #for each gap in the record
    last_flow = df2.loc[df2.index == gaps.start_gap[i], site].item()
    #extract the mean daily flow for the day before the gap starts (.item() ensures this is a number instead    of a dataframe)
    df2.loc[logical_and(df2.index > gaps.start_gap[i], df2.index <= gaps.end_gap[i]), [site]] = last_flow
    #set the flow value during a gap to the value for the final day before the gap starts

df2['Full'] = 0 #create a new column and set all values in it to 0
df2['Full'][df2[site] < fr_flow] = 1 #Set column to 1 where flow is less than full restriction threshold

condition = logical_and(df2[site] < pr_flow, df2[site] >= fr_flow) #find dates where flow is between full and partial restriction threshold

df2['Partial'] = 0
df2['Partial'][condition] = 1

df3 = df2[df2.index.month.isin([9,10,11,12,1,2,3,4])] #new dataframe containing only months in the irrigation season (September - April)

Days_on_restriction = df3.groupby([df3.index.year, df3.index.month])['Full', 'Partial'].sum() #Count the number of days in each month of each year on full and partial restriction

Days_on_restriction1 = Days_on_restriction.unstack() #reformat to table

Days_on_restriction1.rename(columns = {1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'},  inplace=True) #Change column headers to month names

Days_on_restriction1.index.names = ['Year'] #rename row index
Days_on_restriction1.columns.names = ['Restriction Level', 'Month'] #rename column index

export_csv = r'S:\Surface Water\temp\ConwaySH1_Days_on_restriction.csv'
Days_on_restriction1.to_csv(path_or_buf=export_csv)

