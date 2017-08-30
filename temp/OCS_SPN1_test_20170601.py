import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import ephem

s= 'bmh'
plt.style.use(s)

path = '/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/keane/SPN1_20170601/'

### load spn1 data
f1 = 'SPN1_20170524/A1595_Mask_thin.txt'
df1 = pd.read_csv(path+f1,header=4,delimiter='\t',parse_dates=[0],index_col='Time')
df1.Total = pd.to_numeric(df1.Total, errors='coerce')
df1_min = df1.resample('T').mean()

f2 = 'SPN1_20170524/A1596_NOMASK_thin.txt'
df2 = pd.read_csv(path+f2,header=4,delimiter='\t',parse_dates=[0],index_col='Time')
df2.Total = pd.to_numeric(df2.Total, errors='coerce')
df2_min = df2.resample('T').mean()


f3 = 'EppleyData/flex0010-native.swr'
df3 = pd.read_csv(path+f3,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']],index_col='Date_Time')
df3_min = df3.resample('T').mean()

f4 = 'EppleyData/flex1000-native.swr' # truth?
df4 = pd.read_csv(path+f4,header=5,delimiter='\s+', 
		names=['Date','Time', 'Total', 'StdDev', 'Nsec'], 
		na_values = [1E+35], parse_dates=[['Date', 'Time']],index_col='Date_Time')

f5 = 'UW_Radn_rooftop.xlsx'
df5 = pd.read_excel(path+f5,parse_dates=[['Date', 'Time']],index_col='Date_Time')

f6 = 'ISISData/sea17147.dat'
colnames= ['year','jday','month','day','hour','min','dt',
           'zen','dw_psp','qc_dwpsp','direct','qc_direct',
           'diffuse','qc_diffuse','uvb','qc_uvb','uvb_temp',
           'qc_uvb_temp','std_dw_psp','std_direct','std_diffuse',
           'std_uvb']
parse = lambda x: datetime.datetime.strptime(x, '%Y %m %d %H %M')
df6 = pd.read_csv(path+f6,skiprows=2,delimiter='\s+',names=colnames,
                  parse_dates={'datetime':['year','month','day','hour','min']},
                  date_parser=parse,index_col='datetime')

f7 = 'ISISData/sea17148.dat'
colnames= ['year','jday','month','day','hour','min','dt',
           'zen','dw_psp','qc_dwpsp','direct','qc_direct',
           'diffuse','qc_diffuse','uvb','qc_uvb','uvb_temp',
           'qc_uvb_temp','std_dw_psp','std_direct','std_diffuse',
           'std_uvb']
parse = lambda x: datetime.datetime.strptime(x, '%Y %m %d %H %M')
df7 = pd.read_csv(path+f7,skiprows=2,delimiter='\s+',names=colnames,
                  parse_dates={'datetime':['year','month','day','hour','min']},
                  date_parser=parse,index_col='datetime')

plt.figure()
plt.subplot(1,1,1)
df1_min.plot(x=df1_min.index,y=['Total'], ax=plt.gca())
df2_min.plot(x=df2_min.index,y=['Total'], ax=plt.gca())
df2_min.plot(x=df2_min.index,y=['Diffuse'], ax=plt.gca())
df3.plot(x=df3.index,y=['Total'], ax=plt.gca())
df4.plot(x=df4.index,y=['Total'], ax=plt.gca())
df6.plot(x=df6.index,y=['dw_psp'], ax=plt.gca())
df7.plot(x=df7.index,y=['dw_psp'], ax=plt.gca())
ax1=plt.gca()
ax1.legend(['A1595','A1596 Center', 'A1596','flex 0010','flex 1000','UW rooftop'])

df5_s = df5.tshift(-29, freq='T')

###subset
df1_sub = df1_min.loc['2017-05-26 00:00:00':'2017-05-30 00:00:00'].copy()
df2_sub = df2_min.loc['2017-05-26 00:00:00':'2017-05-30 00:00:00'].copy()
df3_sub = df3.loc['2017-05-26 00:00:00':'2017-05-30 00:00:00'].copy()
df4_sub = df4.loc['2017-05-26 00:00:00':'2017-05-30 00:00:00'].copy()
df5_sub = df5.loc['2017-05-26 00:00:00':'2017-05-30 00:00:00']

plt.figure()
plt.subplot(1,1,1)
plt.plot(df4_sub.index,df4_sub['Total']-df1_sub['Total'],'.',markersize=5)
plt.plot(df4_sub.index,df4_sub['Total']-df2_sub['Total'],'.',markersize=5)
plt.plot(df4_sub.index,df4_sub['Total']-df2_sub['Diffuse'],'.',markersize=5)
plt.plot(df4_sub.index,df4_sub['Total']-df3_sub['Total'],'.',markersize=5)
ax1=plt.gca()
ax1.legend(['flex1000-A1595','flex1000-A1596 Center', 'flex1000-A1596','flex1000-flex0010'])
ax1.set_ylim([-200,200])

plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
plt.plot(df4_sub['Total'],df1_sub['Total'],'.')
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


### plot as function of solar angle
def solar_zenith(timestr='2016-01-01 00:00:00', lat='47.6', lon='-122.32'):
    '''
    Time needs to be a string in UTC
    using subroutine ephem for solar positions
    '''
    location = ephem.Observer()
    location.lon, location.lat, location.date = lon,lat,timestr
    sun = ephem.Sun()
    sun.compute(location)
    sun_az = str(sun.az).split(':')
    sun_el = str(sun.alt).split(':')
    sun_az = float(sun_az[0]) + float(sun_az[1])/60. + float(sun_az[2])/3600.
    sun_el = float(sun_el[0]) + float(sun_el[1])/60. + float(sun_el[2])/3600.

    return (90 - sun_el, sun_az) #return solar zenith angle

#seattle
df1_sub['sun_el'] = 0.0
df2_sub['sun_el'] = 0.0
df3_sub['sun_el'] = 0.0
df4_sub['sun_el'] = 0.0
for index, row in df1_sub.iterrows():
    row['sun_el'] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S'),'47.6','-122.32')[0]
for index, row in df2_sub.iterrows():
    row['sun_el'] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S'),'47.6','-122.32')[0]
for index, row in df3_sub.iterrows():
    row['sun_el'] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S'),'47.6','-122.32')[0]
for index, row in df4_sub.iterrows():
    row['sun_el'] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S'),'47.6','-122.32')[0]

plt.figure(figsize=(12,12))
plt.subplot(1,1,1)
plt.plot(df1_sub['sun_el'],df3_sub['Total']/df1_sub['Total'],'.')
plt.plot(df1_sub['sun_el'],df3_sub['Total']/df2_sub['Total'],'.')
plt.plot(df1_sub['sun_el'],df3_sub['Total']/df2_sub['Diffuse'],'.')
plt.plot(df1_sub['sun_el'],df3_sub['Total']/df4_sub['Total'],'.')
plt.ylabel('flex0010/Instrument')
plt.xlabel('Solar Zenith Angle')
plt.ylim([0.75,1.25])
