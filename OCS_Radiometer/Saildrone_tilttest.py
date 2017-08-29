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
				
def rph2za(pitch,roll,heading):
    '''
    from Fortran/IDL subroutine by Paul Ricchiazzi, 25 Feb 97 paul@icess.ucsb.edu

    all angles in degrees
    pitch  : positive values indicate nose up
    roll   : positive values indicate right wing down
    azimuth: positive values clockwise, w.r.t. NORTH
    assume aircraft heading north: rotate around roll
    
    '''
    
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
    '''
    zs = np.cos(np.deg2rad(sunzen))
    ys = np.sin(np.deg2rad(sunzen))*np.cos(np.deg2rad(sunaz))
    xs = np.sin(np.deg2rad(sunzen))*np.sin(np.deg2rad(sunaz))
    zv = np.cos(np.deg2rad(nrmzen))
    yv = np.sin(np.deg2rad(nrmzen))*np.cos(np.deg2rad(nrmaz))
    xv = np.sin(np.deg2rad(nrmzen))*np.sin(np.deg2rad(nrmaz))
    
    muslope = xs*xv + ys*yv + zs*zv
    
    return muslope

"""------------------------------- MAIN ----------------------------------------"""

sdf1='/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/keane/sd_shakedown/sd-1005_20170828T1900-20170829T0400_radiometers-v01.csv'
sd1 = pd.read_csv(sdf1,parse_dates=['isotime'])
sd1.set_index(pd.DatetimeIndex(sd1['isotime']),inplace=True)
sd1_200msec = sd1.resample('200L').mean()
sd1_1sec = sd1.resample('S').mean()

sdf2='/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/keane/sd_shakedown/sd-1006_20170828T1900-20170829T0400_radiometers-v01.csv'
sd2 = pd.read_csv(sdf2,parse_dates=['isotime'])
sd2.set_index(pd.DatetimeIndex(sd2['isotime']),inplace=True)
sd2_200msec = sd2.resample('200L').mean()
sd2_1sec = sd2.resample('S').mean()


"""---"""

"""
fig, axes = plt.subplots(nrows=2, ncols=1)

sd1['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[0]);
sd1['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[0]);
sd1['sw_unshaded_radiometer:raw:average_detector (watt_per_m2)'].plot(ax=axes[0]);
sd1['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'].plot(ax=axes[0]);
axes[0].legend(['sw:raw:shaded total','sw:raw:unshaded center','sw:raw:unshaded ave','sw:raw:shaded diffuse'])
axes[0].set_xticks([])
axes[0].set_title(['sd1005'])


sd3_1sec['lw_radiometer:raw:net_radiation (watt_per_m2)'].plot(ax=axes[1]);
axes[1].legend(['lw:raw:net rad sd1005'])
"""

"""---"""

fig, axes = plt.subplots(nrows=1, ncols=1)
sd1_1sec['lw_radiometer:raw:net_radiation (watt_per_m2)'].plot(ax=axes);
sd2_1sec['lw_radiometer:raw:net_radiation (watt_per_m2)'].plot(ax=axes);
plt.legend(['lw:raw:net rad sd1005','lw:raw:net rad sd1006'])

fig, axes = plt.subplots(nrows=3, ncols=1)
sd1_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[0]);
sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[0]);
axes[0].set_title('sw:raw:total')
axes[0].set_xticks([])
axes[0].legend(['sd1005','sd1006'])

sd1_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[1]);
sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[1]);
axes[1].set_title('sw:raw:unshaded center')
plt.legend(['sd1005','sd1006'])
axes[1].set_xticks([])
axes[1].legend(['sd1005','sd1006'])

sd1_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[2]);
sd2_200msec['sw_shaded_radiometer:raw:total (watt_per_m2)'].plot(ax=axes[2]);
sd1_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[2]);
sd2_200msec['sw_unshaded_radiometer:raw:center_detector (watt_per_m2)'].plot(ax=axes[2]);
sd1_200msec['sw_unshaded_radiometer:raw:average_detector (watt_per_m2)'].plot(ax=axes[2]);
sd2_200msec['sw_unshaded_radiometer:raw:average_detector (watt_per_m2)'].plot(ax=axes[2]);
axes[2].set_title('all uncorrected sw')
axes[2].legend(['sd1005 total','sd1006 total',
                'sd1005 center unshaded','sd1006 center unshaded',
                'sd1005 averaged unshaded','sd1006 averaged unshaded'])


fig, axes = plt.subplots(nrows=3, ncols=1)
sd1_200msec['sw_shaded_radiometer:raw:wing_roll (degrees)'].plot(ax=axes[0]);
sd2_200msec['sw_shaded_radiometer:raw:wing_roll (degrees)'].plot(ax=axes[0]);
axes[0].set_title('sw:raw:roll')
axes[0].set_xticks([])
axes[0].legend(['sd1005','sd1006'])

sd1_200msec['sw_shaded_radiometer:raw:wing_pitch (degrees)'].plot(ax=axes[1]);
sd2_200msec['sw_shaded_radiometer:raw:wing_pitch (degrees)'].plot(ax=axes[1]);
axes[1].set_title('sw:raw:pitch')
axes[1].set_xticks([])
axes[1].legend(['sd1005','sd1006'])

sd1_200msec['sw_shaded_radiometer:raw:wing_yaw (degrees)'].plot(ax=axes[2]);
sd2_200msec['sw_shaded_radiometer:raw:wing_yaw (degrees)'].plot(ax=axes[2]);
axes[2].set_title('sw:raw:yaw')
axes[2].legend(['sd1005','sd1006'])


#### Corrected tilt

sd1_200msec['tilt_normal'] = np.nan
sd1_200msec['sd1_tilt_sza']= np.nan
sd1_200msec['sd1_tilt_az']=  np.nan
sd1_200msec['sd1_tilt_cos_sza']=np.nan

count=0
for index,row in sd1_200msec.iterrows():
    
    if not np.isnan(row['sw_shaded_radiometer:raw:total (watt_per_m2)']) and (count % 5 == 0):
        sd1_200msec['sd1_tilt_sza'][index], sd1_200msec['sd1_tilt_az'][index] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S.%f'),lat=str(row['gps:location:lat (gps_coord)']),lon=str(row['gps:location:lng (gps_coord)']))
  

        instzen,instaz = rph2za(row['sw_shaded_radiometer:raw:pitch (degrees)'],
                                row['sw_shaded_radiometer:raw:roll (degrees)'],
                                row['sw_shaded_radiometer:raw:yaw (degrees)'])
        cos_sza = muslope(sd1_200msec['sd1_tilt_sza'][index],sd1_200msec['sd1_tilt_az'][index],instzen,instaz)
        sd1_200msec['sd1_tilt_cos_sza'][index] = np.rad2deg(np.arccos(cos_sza))
        k_ratio = ( row['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'] * cos_sza ) / (row['sw_shaded_radiometer:raw:total (watt_per_m2)'] - row['sw_shaded_radiometer:raw:diffuse (watt_per_m2)']) 
        sd1_200msec['tilt_normal'][index] = (np.cos(np.deg2rad(sd1_200msec['sd1_tilt_sza'][index])) + k_ratio)/(cos_sza + k_ratio) 
        
    if count % 100 == 0:
        print count
    count +=1  
    
sd2_200msec['tilt_normal'] = np.nan
sd2_200msec['sd2_tilt_sza']= np.nan
sd2_200msec['sd2_tilt_az']=  np.nan
sd2_200msec['sd2_tilt_cos_sza']=np.nan

count=0
for index,row in sd2_200msec.iterrows():
    
    if not np.isnan(row['sw_shaded_radiometer:raw:total (watt_per_m2)']) and (count % 5 == 0):
        sd2_200msec['sd2_tilt_sza'][index], sd2_200msec['sd2_tilt_az'][index] = solar_zenith(index.strftime('%Y-%m-%d %H:%M:%S.%f'),lat=str(row['gps:location:lat (gps_coord)']),lon=str(row['gps:location:lng (gps_coord)']))
  

        instzen,instaz = rph2za(row['sw_shaded_radiometer:raw:pitch (degrees)'],
                                row['sw_shaded_radiometer:raw:roll (degrees)'],
                                row['sw_shaded_radiometer:raw:yaw (degrees)'])
        cos_sza = muslope(sd2_200msec['sd2_tilt_sza'][index],sd2_200msec['sd2_tilt_az'][index],instzen,instaz)
        sd2_200msec['sd2_tilt_cos_sza'][index] = np.rad2deg(np.arccos(cos_sza))
        k_ratio = ( row['sw_shaded_radiometer:raw:diffuse (watt_per_m2)'] * cos_sza ) / (row['sw_shaded_radiometer:raw:total (watt_per_m2)'] - row['sw_shaded_radiometer:raw:diffuse (watt_per_m2)']) 
        sd2_200msec['tilt_normal'][index] = (np.cos(np.deg2rad(sd2['sd2_tilt_sza'][index])) + k_ratio)/(cos_sza + k_ratio) 
        
    if count % 100 == 0:
        print count
    count +=1  