# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 07:35:52 2016

@author: bell
"""

#System Stack
import datetime, sys, os
import argparse

import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='Convert and grid GAK1 historic data')
parser.add_argument('DataPath', metavar='DataPath', type=str,
               help='full path to file')
parser.add_argument('-head','--headerlines',type=int, default=1, 
               help="number of headerlines to skip")
parser.add_argument('-c','--columns',type=str, nargs='+',
               help="space seperated names of colums")
parser.add_argument('-s','--sample',type=str, default='downsample',
               help="upsample/downsample")

args = parser.parse_args()

data = pd.read_csv(args.DataPath,skiprows=args.headerlines,names=args.columns,sep='\s+')

data['DateTime'] =pd.to_datetime(data[['year','month','day','hour','minute']])
data.set_index(pd.DatetimeIndex(data['DateTime']),inplace=True)
data.drop(['year','month','day','hour','minute','second','DateTime'],axis=1,inplace=True)

if args.sample == 'downsample':
    data_hr = data.resample('H', label='right').mean()
elif args.sample == 'upsample':
    data_hr = data.resample('H').interpolate(method='linear')
else:
    data_hr = data

print data_hr.to_csv()