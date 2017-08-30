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

path = '/Volumes/WDC_internal/Users/bell/Data_Local/SPN1_Development/OCS_SPN1Tests/ArcticSensorComparison/'

### load spn1 data
f1 = 'SPN1_A1274_Lakeside.txt'
df1 = pd.read_csv(path+f1,header=4,delimiter='\t')
df1['Time'] = pd.to_datetime(df1['Time'])
df1['Total'] = pd.to_numeric(df1['Total'], errors='coerce')
df1['Direct'] =df1['Total']-df1['Diffuse']

f2 = 'SPN1_A1596_Lakeside.txt'
df2 = pd.read_csv(path+f2,header=4,delimiter='\t')
df2['Time'] = pd.to_datetime(df2['Time'])
df2['Direct'] =df2['Total']-df2['Diffuse']

###plot 4 instruments
# SPN1, SPN1, Eppley

plt.figure()
plt.subplot(2,1,1)
df2.plot(x='Time',y=['Total'], ax=plt.gca())
df1.plot(x='Time',y=['Total'], ax=plt.gca())
plt.legend([f2,f1],loc='upper center',ncol=4)
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=[0,6,12,18]))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)

plt.subplot(2,1,2)
df2.plot(x='Time',y=['Diffuse'], ax=plt.gca())
df1.plot(x='Time',y=['Diffuse'], ax=plt.gca())
plt.legend([f2,f1])
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=[0,6,12,18]))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)


### gridded
# map all variables to same time
ts = datetime.datetime(2016,11,8)
numdays = 24*60*2
date_list = [ts + datetime.timedelta(0,0,0,0,x) for x in range(0, numdays)]
SPN1_A1274d = collections.OrderedDict()
SPN1_A1274difd = collections.OrderedDict()
SPN1_A1596d= collections.OrderedDict()
SPN1_A1596difd = collections.OrderedDict()
for date_val in date_list:
	if (date_val < datetime.datetime(2016,8,13)) or (date_val >= datetime.datetime(2016,8,14)):
		if (df1['Time']==date_val).any():
			df1_ind = (df1['Time']==date_val).argmax()
			SPN1_A1274 = df1.iloc[df1_ind]['Total']
			SPN1_A1274dif = df1.iloc[df1_ind]['Diffuse']
		else:
			SPN1_A1274 = np.nan
			SPN1_A1274dif = np.nan

		if (df2['Time']==date_val).any():
			df2_ind = (df2['Time']==date_val).argmax()
			SPN1_A1596 = df2.iloc[df2_ind]['Total']
			SPN1_A1596dif = df2.iloc[df2_ind]['Diffuse']
		else:
			SPN1_A1596 = np.nan
			SPN1_A1596dif = np.nan



		SPN1_A1274d[date_val] = SPN1_A1274
		SPN1_A1274difd[date_val] = SPN1_A1274dif
		SPN1_A1596d[date_val] = SPN1_A1596
		SPN1_A1596difd[date_val] = SPN1_A1596dif


### Calculate Direct from SPN1
plt.figure()
plt.subplot(2,1,1)
df1.plot(ax=plt.gca())
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=[0,6,12,18]))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)
plt.subplot(2,1,2)
df2.plot(ax=plt.gca())
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=[0,6,12,18]))
ax1.xaxis.set_minor_locator(HourLocator(byhour=12))
ax1.xaxis.set_major_formatter(ticker.NullFormatter())
ax1.xaxis.set_minor_formatter(DateFormatter('%b %d'))
ax1.tick_params(axis='both', which='minor', labelsize=12)

###
plt.figure()
plt.subplot(2,1,1)
(df1[["P1","P2","P3","P4","P5","P6","P7"]].min(axis=1)*2.0*1.02*1.14).plot(ax=plt.gca())
(df1[["P1","P2","P3","P4","P5","P6","P7"]].max(axis=1)*2.0*1.02*1.14).plot(ax=plt.gca())
(df1[["P1","P2","P3","P4","P5","P6","P7"]]*2.0*1.02*1.14).std(axis=1).plot(ax=plt.gca())
df1[["Diffuse"]].plot(ax=plt.gca())
df1[["Total"]].plot(ax=plt.gca())
plt.legend(['7ch min diffuse','7ch max diffuse','7ch std','Diffuse (spn1 calc)','Total (spn1 calc)'])
plt.subplot(2,1,2)
(df2[["P1","P2","P3","P4","P5","P6","P7"]].min(axis=1)*2.0*1.02*1.14).plot(ax=plt.gca())
(df2[["P1","P2","P3","P4","P5","P6","P7"]].max(axis=1)*2.0*1.02*1.14).plot(ax=plt.gca())
(df2[["P1","P2","P3","P4","P5","P6","P7"]]*2.0*1.02*1.14).std(axis=1).plot(ax=plt.gca())
df2[["Diffuse"]].plot(ax=plt.gca())
df2[["Total"]].plot(ax=plt.gca())
plt.legend(['7ch min diffuse','7ch max diffuse','7ch std','Diffuse (spn1 calc)','Total (spn1 calc)'])
