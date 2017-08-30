# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 08:13:16 2016

@author: bell
"""

import os
import pandas as pd

path = '/Users/bell/Data_Local/SPN1_Development/ukmet_buoy/download/L4/SPN1/'
all_files= [each for each in os.listdir(path) if each.endswith('.dat')]

frame = pd.DataFrame()

for f in all_files[610:620]:
    df = pd.read_csv(path+f,header=0,names=['date', 'time', 'total', 'diffuse', 'sunshinebit'])
    df['total'] = df['total'].str[2:].astype(float)
    df['date'] = pd.to_datetime(df['date']+df['time'], dayfirst=True)
    frame = frame.append(df, ignore_index=True)

frame.set_index('date',inplace=True)

frame.plot()
#pd.rolling_std(frame['total'], 10).plot()
#pd.rolling_mean(frame['total'],6).plot()


"""
 spn_file = '/Users/bell/Programs/Python/MooringDataProcessing/oer_radiometer/data/deploy_2015/2015_OER_SWR_SPN1_0904.csv'
df2 = pd.read_csv(spn_file)
df2['SPN Total 1min'] = pd.rolling_mean(df2['SPN Total'],12)
df2['SPN Diffuse 1min'] = pd.rolling_mean(df2['SPN Diffuse'],12)
df2.ix[:, ['SPN Total 1min','SPN Diffuse 1min','SWR']].plot()
df2.ix[:, ['SPN Total','SPN Diffuse','SWR']].plot()

pd.rolling_std(df2['SPN Total'], 10).plot()
"""