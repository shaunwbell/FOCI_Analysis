#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:03:49 2017

@author: bell
"""

import xarray as xa
import pandas as pd
import datetime
from netCDF4 import num2date

"""
filein='/Volumes/WDC_internal/Users/bell/in_and_outbox/Ongoing_Analysis/M2_IntegratedTemp/most_current/1995_2016_htcontent_depthint.cf.nc'
data = xa.open_dataset(filein)
df = data.to_dataframe()
df=df.reset_index().set_index('time')
for month in range(1,13,1):
	for day in range(1,32):
		temp = [str(month)+'-'+str(day)]
		for year in range(1995,2017):
			temp = temp + [df.loc[(df.index.day == day) & (df.index.month == month) & (df.index.year == year)].mean()['V00_1900']]
		print temp
        
"""        

filein='/Volumes/WDC_internal/Users/bell/ecoraid/2016/Moorings/16bspitae/initial_archive/16bspitae_prawler_2D.nc'
data = xa.open_dataset(filein, decode_cf=False)
datamasked = data.where(data.time < 1e30)

time = datamasked.time.to_pandas()
time = time.fillna(0)
tarray = num2date(time,'hours since 1900-01-01T00:00:00Z')

for i in range(0, tarray.shape[0]):
    line = []
    for j in range(0, tarray.shape[1]):
        if tarray[i,j] != datetime.datetime(1900,1,1):
            line = line + [datetime.datetime.strftime(tarray[i,j],'%Y-%m-%d %H:%M:%S')]
        else:
            line = line + ['']
    print ",".join(line)
#print Chlorophyll