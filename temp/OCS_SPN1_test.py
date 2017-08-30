# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 07:35:52 2016

@author: bell
"""

import pandas as pd
import numpy as np
import collections
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, WeekdayLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import matplotlib.ticker as ticker

s= 'bmh'
plt.style.use(s)

path = '/Users/bell/Downloads/SPN1_Development/OCS_SPN1Tests/GK_1minAvg/'
pathb = '/Users/bell/Downloads/SPN1_Development/OCS_SPN1Tests/'

### load spn1 data
f1 = 'SPN1_A1593.txt'
df1 = pd.read_csv(path+f1,header=4,delimiter='\t')
df1['Time'] = pd.to_datetime(df1['Time'])

f2 = 'SPN1_A1594.txt'
df2 = pd.read_csv(path+f2,header=4,delimiter='\t')
df2['Time'] = pd.to_datetime(df2['Time'])


### load spn1 data 5s
f1b = 'SPN1_A1593_5s.txt'
df1b = pd.read_csv(pathb+f1b,header=3,delimiter='\t')
df1b['Time'] = pd.to_datetime(df1b['Time'])

f2b = 'SPN1_A1594_5s.txt'
df2b = pd.read_csv(pathb+f2b,header=3,delimiter='\t')
df2b['Time'] = pd.to_datetime(df2b['Time'])

f3 = 'flex1000-native.swr'
df3 = pd.read_csv(path+f3,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']])


f4 = 'flex1011-native.swr'
df4 = pd.read_csv(path+f4,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']])

###plot 4 instruments
# SPN1, SPN1, Eppley
plt.figure()
plt.subplot(2,1,1)
df2b.plot(x='Time',y=['Total'], ax=plt.gca())
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=0))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)
plt.subplot(2,1,2)
df1b.plot(x='Time',y=['Diffuse'], ax=plt.gca())
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=0))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)


plt.figure()
plt.subplot(2,1,1)
df2.plot(x='Time',y=['Total'], ax=plt.gca())
df1.plot(x='Time',y=['Total'], ax=plt.gca())
df3.plot(x='Date_Time',y=['Total'], ax=plt.gca())
df4.plot(x='Date_Time',y=['Total'], ax=plt.gca())
plt.legend([f2,f1,f3,f4],loc='upper center',ncol=4)
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=0))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(2,1,2)
df2.plot(x='Time',y=['Diffuse'], ax=plt.gca())
df1.plot(x='Time',y=['Diffuse'], ax=plt.gca())
plt.legend([f2,f1])
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=0))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)

### gridded
# map all variables to same time
ts = datetime.datetime(2016,8,12)
numdays = 24*60*8
date_list = [ts + datetime.timedelta(0,0,0,0,x) for x in range(0, numdays)]
SPN1_A1593d = collections.OrderedDict()
SPN1_A1593difd = collections.OrderedDict()
SPN1_A1594d= collections.OrderedDict()
SPN1_A1594difd = collections.OrderedDict()
flex1000d = collections.OrderedDict()
flex1011d = collections.OrderedDict()
for date_val in date_list:
	if (date_val < datetime.datetime(2016,8,13)) or (date_val >= datetime.datetime(2016,8,14)):
		if (df1['Time']==date_val).any():
			df1_ind = (df1['Time']==date_val).argmax()
			SPN1_A1593 = df1.iloc[df1_ind]['Total']
			SPN1_A1593dif = df1.iloc[df1_ind]['Diffuse']
		else:
			SPN1_A1593 = np.nan
			SPN1_A1593dif = np.nan

		if (df2['Time']==date_val).any():
			df2_ind = (df2['Time']==date_val).argmax()
			SPN1_A1594 = df2.iloc[df2_ind]['Total']
			SPN1_A1594dif = df2.iloc[df2_ind]['Diffuse']
		else:
			SPN1_A1594 = np.nan
			SPN1_A1594dif = np.nan

		if (df3['Date_Time']==date_val).any():
			df3_ind = (df3['Date_Time']==date_val).argmax()
			flex1000 = df3.iloc[df3_ind]['Total']
		else:
			flex1000 = np.nan

		if (df4['Date_Time']==date_val).any():
			df4_ind = (df4['Date_Time']==date_val).argmax()
			flex1011 = df4.iloc[df4_ind]['Total']
		else:
			flex1011 = np.nan

		SPN1_A1593d[date_val] = SPN1_A1593
		SPN1_A1593difd[date_val] = SPN1_A1593dif
		SPN1_A1594d[date_val] = SPN1_A1594
		SPN1_A1594difd[date_val] = SPN1_A1594dif
		flex1000d[date_val] = flex1000
		flex1011d[date_val] = flex1011

### Calculate Direct from SPN1
plt.figure()
plt.subplot(2,1,1)
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593d.values()),'.k')
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593d.values()) - np.array(SPN1_A1593difd.values()),'.r')
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593difd.values()),'.b')
plt.subplot(2,1,2)
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594d.values()),'.k')
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594d.values()) - np.array(SPN1_A1594difd.values()),'.r')
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594difd.values()),'.b')

### Difference measurements
plt.figure()
plt.subplot(4,1,1)
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593d.values()) - np.array(SPN1_A1594d.values()),'.k')
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593d.values()) - np.array(flex1000d.values()),'.r')
plt.plot(SPN1_A1593d.keys(), np.array(SPN1_A1593d.values()) - np.array(flex1011d.values()),'.b')
plt.subplot(4,1,2)
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594d.values()) - np.array(SPN1_A1593d.values()),'.g')
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594d.values()) - np.array(flex1000d.values()),'.r')
plt.plot(SPN1_A1594d.keys(), np.array(SPN1_A1594d.values()) - np.array(flex1011d.values()),'.b')
plt.subplot(4,1,3)
plt.plot(flex1000d.keys(), np.array(flex1000d.values()) - np.array(SPN1_A1593d.values()),'.g')
plt.plot(flex1000d.keys(), np.array(flex1000d.values()) - np.array(SPN1_A1594d.values()),'.k')
plt.plot(flex1000d.keys(), np.array(flex1000d.values()) - np.array(flex1011d.values()),'.b')
plt.subplot(4,1,4)
plt.plot(flex1011d.keys(), np.array(flex1011d.values()) - np.array(SPN1_A1593d.values()),'.g')
plt.plot(flex1011d.keys(), np.array(flex1011d.values()) - np.array(SPN1_A1594d.values()),'.k')
plt.plot(flex1011d.keys(), np.array(flex1011d.values()) - np.array(flex1000d.values()),'.r')

### Grid comparison of each instrument
plt.figure()
plt.subplot(4,4,1)
plt.plot(SPN1_A1593d.values(),SPN1_A1593d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,2)
plt.plot(SPN1_A1594d.values(),SPN1_A1593d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,3)
plt.plot(flex1000d.values(),SPN1_A1593d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,4)
plt.plot(flex1011d.values(),SPN1_A1593d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,5)
plt.plot(SPN1_A1593d.values(),SPN1_A1594d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,6)
plt.plot(SPN1_A1594d.values(),SPN1_A1594d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,7)
plt.plot(flex1000d.values(),SPN1_A1594d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,8)
plt.plot(flex1011d.values(),SPN1_A1594d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,9)
plt.plot(SPN1_A1593d.values(),flex1000d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,10)
plt.plot(SPN1_A1594d.values(),flex1000d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,11)
plt.plot(flex1000d.values(),flex1000d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,12)
plt.plot(flex1011d.values(),flex1000d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,13)
plt.plot(SPN1_A1593d.values(),flex1011d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,14)
plt.plot(SPN1_A1594d.values(),flex1011d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,15)
plt.plot(flex1000d.values(),flex1011d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')

plt.subplot(4,4,16)
plt.plot(flex1011d.values(),flex1011d.values(),'.k')
plt.plot(range(0,1000,100),range(0,1000,100),'-r')
