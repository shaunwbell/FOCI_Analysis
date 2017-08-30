import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, WeekdayLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import matplotlib.ticker as ticker

s= 'bmh'
plt.style.use(s)

path = '/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/keane/SPN1_20170518/'

### load spn1 data
f1 = 'A1595_Mask.txt'
df1 = pd.read_csv(path+f1,header=4,delimiter='\t',parse_dates=[0])
df1 = df1.set_index(pd.DatetimeIndex(df1['Time']))
df1.Total = pd.to_numeric(df1.Total, errors='coerce')
df1_min = df1.resample('T').mean()

f2 = 'A1596_NOMASK.txt'
df2 = pd.read_csv(path+f2,header=3,delimiter='\t',parse_dates=[0])
df2 = df2.set_index(pd.DatetimeIndex(df2['Time']))
df2_min = df2.resample('T').mean()

f3 = 'flex0010-native.swr'
df3 = pd.read_csv(path+f3,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']])
df3 = df3.set_index(pd.DatetimeIndex(df3['Date_Time']))

f4 = 'flex1000-native.swr' # truth?
df4 = pd.read_csv(path+f4,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']])
df4 = df4.set_index(pd.DatetimeIndex(df4['Date_Time']))

f5 = 'UW_Radn_rooftop.xlsx'
df5 = pd.read_excel(path+f5,parse_dates=[['Date', 'Time']])

plt.figure()
plt.subplot(1,1,1)
df1_min.plot(x=df1_min.index,y=['Total'], ax=plt.gca())
df2_min.plot(x=df2_min.index,y=['Total'], ax=plt.gca())
df2_min.plot(x=df2_min.index,y=['Diffuse'], ax=plt.gca())
df3.plot(x=df3.index,y=['Total'], ax=plt.gca())
df4.plot(x=df4.index,y=['Total'], ax=plt.gca())
df5.plot(x='Date_Time',y=['Radn'], ax=plt.gca())
ax1=plt.gca()
ax1.legend(['A1595','A1596 Center', 'A1596','flex 0010','flex 1000','UW rooftop'])

###subset
df1_sub = df1_min.loc['2017-05-20 20:00:00':'2017-05-21 17:00:00']
df2_sub = df2_min.loc['2017-05-20 20:00:00':'2017-05-21 17:00:00']
df3_sub = df3.loc['2017-05-20 20:00:00':'2017-05-21 17:00:00']
df4_sub = df4.loc['2017-05-20 20:00:00':'2017-05-21 17:00:00']

plt.figure()
plt.subplot(1,1,1)
plt.plot(df4_sub.index,df4_sub['Total']-df1_sub['Total'],'.')
plt.plot(df4_sub.index,df4_sub['Total']-df2_sub['Total'],'.')
plt.plot(df4_sub.index,df4_sub['Total']-df2_sub['Diffuse'],'.')
plt.plot(df4_sub.index,df4_sub['Total']-df3_sub['Total'],'.')
ax1=plt.gca()
ax1.legend(['flex1000-A1595','flex1000-A1596 Center', 'flex1000-A1596','flex1000-flex0010'])
ax1=plt.gca()
ax1.xaxis.set_major_locator(HourLocator(byhour=0))
ax1.xaxis.set_minor_locator(HourLocator(byhour=range(3,24,3)))
ax1.xaxis.set_major_formatter(DateFormatter('%b %d'))
ax1.xaxis.set_minor_formatter(DateFormatter('%H:%M:%S'))
ax1.tick_params(axis='both', which='both', labelsize=12)

plt.figure()
plt.subplot(2,2,1)
plt.plot(df2_sub['Diffuse'],df2_sub['Total'],'.')
plt.plot(range(0,1000,10),range(0,1000,10),'--r')
plt.ylabel('A1595')
plt.xlabel('flex1000')
plt.subplot(2,2,2)
plt.plot(df4_sub['Total'],df2_sub['Total'],'.')
plt.plot(range(0,1000,10),range(0,1000,10),'--r')
plt.ylabel('A1596 Center')
plt.xlabel('flex1000')
plt.subplot(2,2,3)
plt.plot(df4_sub['Total'],df2_sub['Diffuse'],'.')
plt.plot(range(0,1000,10),range(0,1000,10),'--r')
plt.ylabel('A1596')
plt.xlabel('flex1000')
plt.subplot(2,2,4)
plt.plot(df4_sub['Total'],df3_sub['Total'],'.')
plt.plot(range(0,1000,10),range(0,1000,10),'--r')
plt.ylabel('flex0010')
plt.xlabel('flex1000')
