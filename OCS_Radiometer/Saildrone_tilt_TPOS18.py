#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 07:01:14 2017

@author: bell
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import ephem
def solar_zenith(timestr='2016-01-01 00:00:00', lat='47.6', lon='-122.32'):
    '''
    Time needs to be a string in UTC
    using subroutine ephem for solar positions
    '''
    location = ephem.Observer()
    location.lon, location.lat, location.date = lon,lat,timestr
    sun = ephem.Sun()
    sun.compute(location)
    
    sun_string = True
    if sun_string: 
        sun_az = str(sun.az).split(':')
        sun_el = str(sun.alt).split(':')
        sun_az = float(sun_az[0]) + float(sun_az[1])/60. + float(sun_az[2])/3600.
        sun_el = float(sun_el[0]) + float(sun_el[1])/60. + float(sun_el[2])/3600.
    else:
        sun_az = ephem.degrees(sun.az)
        sun_el = ephem.degrees(sun.alt)
        
    return (90. - sun_el, sun_az)
				
def prh2za(pitch,roll,heading):
    '''
    from Fortran/IDL subroutine by Paul Ricchiazzi, 25 Feb 97 paul@icess.ucsb.edu

    all angles in degrees
    pitch  : positive values indicate nose up
    roll   : positive values indicate right wing down
    azimuth: positive values clockwise, w.r.t. NORTH
    assume aircraft heading north: rotate around roll

    Input: Either Pandas Series or Numpy Array
    
    '''
    
    if isinstance(pitch, pd.Series):
        ### pandas series
        roll_rad = np.deg2rad(roll)
        pitch_rad = np.deg2rad(pitch)
        heading_rad = np.deg2rad(heading)

        uz = np.cos(roll_rad)*np.cos(pitch_rad)
        ux = np.sin(roll_rad)
        uy = -1.0 * np.cos(roll_rad) * np.sin(pitch_rad)

        vz = uz
        vx = ux * np.cos(heading_rad) + uy * np.sin(heading_rad)
        vy = uy * np.cos(heading_rad) - uy * np.sin(heading_rad)

        zenith = np.rad2deg(np.arccos(vz))
        vy.loc[ vy == 0.0 ] = 0.00000001
        azimuth = np.rad2deg(np.arctan2(vx,vy))
        
        return (zenith.values, azimuth.values)

    else:
        roll_rad = np.deg2rad(np.float(roll))
        pitch_rad = np.deg2rad(np.float(pitch))
        heading_rad = np.deg2rad(np.float(heading))

        uz = np.cos(roll_rad)*np.cos(pitch_rad)
        ux = np.sin(roll_rad)
        uy = -1.0 * np.cos(roll_rad) * np.sin(pitch_rad)
        
        vz = uz
        vx = ux * np.cos(heading_rad) + uy * np.sin(heading_rad)
        vy = uy * np.cos(heading_rad) - uy * np.sin(heading_rad)
        
        zenith = np.rad2deg(np.arccos(vz))
        if (vy == 0.0):
            vy = 0.00000001
        azimuth = np.rad2deg(np.arctan2(vx,vy))
    
    return (zenith, azimuth)
    
def muslope(sunzen,sunaz,nrmzen,nrmaz):
    '''
    C ROUTINE:  muslope
    C
    C PURPOSE:  compute dot product of surface normal to incomming solar ray
    C
    C USEAGE:   result=muslope(sunzen,sunaz,nrmzen,nrmaz)
    C
    C INPUT:    
    C   sunzen  solar zenith angle (degrees)
    C
    C   sunaz   solar azimuth angle (clockwise from due north, degrees) 
    C
    C   nrmzen  zenith angle of surface normal vector
    C           (nrmzen=0 for a horizontal surface)
    C
    C   nrmaz   azimuth angle of surface normal vector (nrmaz=45 degrees
    C           for a surface that slopes down in the north-east direction)
    C
    C OUTPUT:
    C   result  cosine of angle between sun and instrument 
    C           = cos(sza) if instrument pointing to zenith
    C
    C AUTHOR:   Yan Shi, PNNL, 10/2/2008 
    C           converted from IDL code by Paul Ricchiazzi
    C           Institute for Computational Earth System Science
    C           University of California, Santa Barbara
    C           paul@icess.ucsb.edu

    Input: Either Pandas Series or Numpy Array

    '''


    zs = np.cos(np.deg2rad(sunzen))
    ys = np.sin(np.deg2rad(sunzen))*np.cos(np.deg2rad(sunaz))
    xs = np.sin(np.deg2rad(sunzen))*np.sin(np.deg2rad(sunaz))
    zv = np.cos(np.deg2rad(nrmzen))
    yv = np.sin(np.deg2rad(nrmzen))*np.cos(np.deg2rad(nrmaz))
    xv = np.sin(np.deg2rad(nrmzen))*np.sin(np.deg2rad(nrmaz))
    
    muslope = xs*xv + ys*yv + zs*zv
 
    if isinstance(sunzen, pd.Series):
        ### pandas series    
        return muslope.values
    else:
        return muslope

"""------------------------------- MAIN ----------------------------------------"""

#%%
# Read Data Files
# sd1 is 10s 1005 data
# sd2 is 200ms 1005 data
# sd3 is 200ms 1006 data

sdf1='/Users/bell/in_and_outbox/2018/saildrone_TPOS/sd-1005_20171001T0600-20171002T0600_sw_radiometer-v01.csv'
sd1 = pd.read_csv(sdf1,parse_dates=['isotime'])
sd1.set_index(pd.DatetimeIndex(sd1['isotime']),inplace=True)
#sd1_200msec = sd1.resample('200L').mean()
sd1_1min = sd1.resample('60S').mean()
sd1_10sec = sd1.resample('10S').mean()

sdf2='/Users/bell/in_and_outbox/2018/saildrone_TPOS/sd-1005_20171001T0930-20171001T2030_sw_radiometer-v01.csv'
sd2 = pd.read_csv(sdf2,parse_dates=['isotime'])
sd2.set_index(pd.DatetimeIndex(sd2['isotime']),inplace=True)
sd2_200msec = sd1.resample('200L').mean()
sd2_10sec = sd2.resample('10S').mean()

sd1_roll_offset = -1.0
sd1_pitch_offset = 8.3

sd2_roll_offset = -1.0
sd2_pitch_offset = 8.3

#sd2_roll_offset = -1.2
#sd2_pitch_offset = 7.8
#%%
"""

Saildrone uses NED coordinates...

North east down convention

So positive pitch is nose down, positive yaw is clockwise, positive roll is right wing down

"""

##Longwave
fig, axes = plt.subplots(nrows=1, ncols=1)
sd1_10sec['lw_radiometer:raw:net_radiation (watt_per_m2)'].plot(ax=axes);
sd2_200msec['lw_radiometer:raw:net_radiation (watt_per_m2)'].plot(ax=axes);
sd1_10sec['lw_radiometer:raw:downwelling_radiation (watt_per_m2)'].plot(ax=axes);
sd2_200msec['lw_radiometer:raw:downwelling_radiation (watt_per_m2)'].plot(ax=axes);
plt.legend(['lw:raw:net rad sd1005','lw:raw:net rad sd1006','lw:raw:downwelling rad sd1005','lw:raw:downwelling rad sd1006'])

#%%

##Shortwave
fig, axes = plt.subplots(nrows=3, ncols=1)
sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[0]);
sd2_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[0]);
axes[0].set_title('sw:raw:total')
axes[0].set_xticks([])
axes[0].legend(['sd1005','sd1006'])

sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[1]);
sd2_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[1]);
axes[1].set_title('sw:raw:unshaded center')
plt.legend(['sd1005','sd1006'])
axes[1].set_xticks([])
axes[1].legend(['sd1005','sd1006'])

sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[2]);
sd2_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[2]);
sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[2]);
sd2_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[2]);
sd1_10sec['sw_unshaded_radiometer:raw:average_detector (watt_per_m2)'].plot(ax=axes[2]);
sd2_10sec['sw_unshaded_radiometer:raw:average_detector (watt_per_m2)'].plot(ax=axes[2]);
axes[2].set_title('all uncorrected sw')
axes[2].legend(['sd1005 total','sd1006 total',
                'sd1005 center unshaded','sd1006 center unshaded',
                'sd1005 averaged unshaded','sd1006 averaged unshaded'])

#%%
##platform tilt/roll
fig, axes = plt.subplots(nrows=3, ncols=1)
sd1_10sec['sw_shaded_radiometer:raw:wing_roll (degrees)'].plot(ax=axes[0]);
sd2_10sec['sw_shaded_radiometer:raw:wing_roll (degrees)'].plot(ax=axes[0]);
axes[0].set_title('sw:raw:roll')
axes[0].set_xticks([])
axes[0].legend(['sd1005','sd1006'])

(sd1_10sec['sw_shaded_radiometer:raw:wing_pitch (degrees)']).plot(ax=axes[1]);
(sd2_10sec['sw_shaded_radiometer:raw:wing_pitch (degrees)']).plot(ax=axes[1]);
axes[1].set_title('sw:raw:pitch')
axes[1].set_xticks([])
axes[1].legend(['sd1005','sd1006'])

sd1_10sec['sw_shaded_radiometer:raw:wing_yaw (degrees)'].plot(ax=axes[2]);
sd2_10sec['sw_shaded_radiometer:raw:wing_yaw (degrees)'].plot(ax=axes[2]);
axes[2].set_title('sw:raw:yaw')
axes[2].legend(['sd1005','sd1006'])

#%%
#### Corrected tilt -10s

sd1_10sec['sd1_tilt_az'] = np.nan
sd1_10sec['sd1_tilt_sza']=  np.nan

count=0
for index,row in sd1_10sec.iterrows():
    
    if not np.isnan(row['sw_shaded_radiometer:raw:total (watt_per_m2)']):
        sd1_10sec['sd1_tilt_sza'][index], sd1_10sec['sd1_tilt_az'][index] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S.%f'),lat=str(row['sw_shaded_radiometer:raw:gps_lat (degrees)']),lon=str(row['sw_shaded_radiometer:raw:gps_lng (degrees)']))
    if count % 100 == 0:
        print count
    count +=1   

sd1_10sec['instzen'],sd1_10sec['instaz'] = prh2za(-1*(sd1_pitch_offset+sd1_10sec['sw_shaded_radiometer:raw:wing_pitch (degrees)']),
                        sd1_roll_offset+sd1_10sec['sw_shaded_radiometer:raw:wing_roll (degrees)'],
                        sd1_10sec['sw_shaded_radiometer:raw:wing_yaw (degrees)'])

sd1_10sec['cos_sza'] = muslope(sd1_10sec['sd1_tilt_sza'],sd1_10sec['sd1_tilt_az'],sd1_10sec['instzen'],sd1_10sec['instaz'])
sd1_10sec['sd1_tilt_cos_sza'] = np.rad2deg(np.arccos(sd1_10sec['cos_sza']))
sd1_10sec['k_ratio'] = ( sd1_10sec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'] * sd1_10sec['cos_sza'] ) / (sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'] - sd1_10sec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)']) 
sd1_10sec['tilt_normal'] = (np.cos(np.deg2rad(sd1_10sec['sd1_tilt_sza'])) + sd1_10sec['k_ratio'])/(sd1_10sec['cos_sza'] + sd1_10sec['k_ratio']) 
        
#%%
## Corrected Tilt Raw
sd2_200msec['sd2_tilt_az'] = np.nan
sd2_200msec['sd2_tilt_sza']=  np.nan

count=0
for index,row in sd2_200msec.iterrows():
    
    if not np.isnan(row['sw_shaded_radiometer:raw:total (watt_per_m2)']):
        sd2_200msec['sd2_tilt_sza'][index], sd2_200msec['sd2_tilt_az'][index] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S.%f'),lat=str(row['sw_shaded_radiometer:raw:gps_lat (degrees)']),lon=str(row['sw_shaded_radiometer:raw:gps_lng (degrees)']))
    if count % 100 == 0:
        print count
    count +=1   

sd2_200msec['instzen'],sd2_200msec['instaz'] = prh2za(-1*(sd2_pitch_offset+sd2_200msec['sw_shaded_radiometer:raw:wing_pitch (degrees)']),
                        sd2_roll_offset+sd2_200msec['sw_shaded_radiometer:raw:wing_roll (degrees)'],
                        sd2_200msec['sw_shaded_radiometer:raw:wing_yaw (degrees)'])

sd2_200msec['cos_sza'] = muslope(sd2_200msec['sd2_tilt_sza'],sd2_200msec['sd2_tilt_az'],sd2_200msec['instzen'],sd2_200msec['instaz'])
sd2_200msec['sd2_tilt_cos_sza'] = np.rad2deg(np.arccos(sd2_200msec['cos_sza']))
sd2_200msec['k_ratio'] = ( sd2_200msec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'] * sd2_200msec['cos_sza'] ) / (sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'] - sd2_200msec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)']) 
sd2_200msec['tilt_normal'] = (np.cos(np.deg2rad(sd2_200msec['sd2_tilt_sza'])) + sd2_200msec['k_ratio'])/(sd2_200msec['cos_sza'] + sd2_200msec['k_ratio']) 


#%%        
###plots
   
fig, axes = plt.subplots(nrows=4, ncols=1)
(sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[0])
(sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[0])
(sd1_10sec['tilt_normal'].resample('1s').mean()*sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[0])
(sd1_10sec['tilt_normal'].resample('1s').mean()*sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[0])
axes[0].set_title('sd 1005 tilt corrected')
axes[0].set_xticks([])
axes[0].set_ylim([-100, 1000])
axes[0].legend(['unshaded uncorrected','shaded uncorrected', 'unshaded corrected','shaded corrected'])

(sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[1])
(sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[1])
(sd2_200msec['tilt_normal'].resample('1s').mean()*sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[1])
(sd2_200msec['tilt_normal'].resample('1s').mean()*sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[1])
axes[1].set_title('sd 1006 tilt corrected')
axes[1].set_xticks([])
axes[1].set_ylim([-100, 1000])
axes[1].legend(['unshaded uncorrected','shaded uncorrected', 'unshaded corrected','shaded corrected'])

(sd1_10sec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[2])
(sd2_200msec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[2])
axes[2].set_title('diffuse')
axes[2].set_xticks([])
axes[2].legend(['sd1005','sd1006'])

(sd1_10sec['sd1_tilt_cos_sza'].resample('1s').mean()).plot(ax=axes[3])
(sd2_200msec['sd2_tilt_cos_sza'].resample('1s').mean()).plot(ax=axes[3])
(sd1_10sec['sd1_tilt_sza'].resample('1s').mean()).plot(ax=axes[3])
axes[3].set_title('solar angles')
axes[3].legend(['sd1005 sza','sd1006 sza','local sza'])


fig, axes = plt.subplots(nrows=4, ncols=1)
(sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('10s').mean()).plot(ax=axes[0])
(sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('10s').mean()).plot(ax=axes[0])
(sd1_10sec['tilt_normal'].resample('1s').mean()*sd1_10sec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).resample('10s').mean().plot(ax=axes[0])
(sd1_10sec['tilt_normal'].resample('1s').mean()*sd1_10sec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).resample('10s').mean().plot(ax=axes[0])
axes[0].set_title('sd 1005 tilt corrected')
axes[0].set_xticks([])
axes[0].set_ylim([-100, 1000])
axes[0].legend(['unshaded uncorrected','shaded uncorrected', 'unshaded corrected','shaded corrected'])

(sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('10s').mean()).plot(ax=axes[1])
(sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('10s').mean()).plot(ax=axes[1])
(sd2_200msec['tilt_normal'].resample('1s').mean()*sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].resample('1s').mean()).resample('10s').mean().plot(ax=axes[1])
(sd2_200msec['tilt_normal'].resample('1s').mean()*sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].resample('1s').mean()).resample('10s').mean().plot(ax=axes[1])
axes[1].set_title('sd 1006 tilt corrected')
axes[1].set_xticks([])
axes[1].set_ylim([-100, 1000])
axes[1].legend(['unshaded uncorrected','shaded uncorrected', 'unshaded corrected','shaded corrected'])

(sd1_10sec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[2])
(sd2_200msec['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'].resample('1s').mean()).plot(ax=axes[2])
axes[2].set_title('diffuse')
axes[2].set_xticks([])
axes[2].legend(['sd1005','sd1006'])

(sd1_10sec['sd1_tilt_cos_sza'].resample('1s').mean()).plot(ax=axes[3])
(sd2_200msec['sd2_tilt_cos_sza'].resample('1s').mean()).plot(ax=axes[3])
(sd1_10sec['sd1_tilt_sza'].resample('1s').mean()).plot(ax=axes[3])
axes[3].set_title('solar angles')
axes[3].legend(['sd1005 sza','sd1006 sza','local sza'])