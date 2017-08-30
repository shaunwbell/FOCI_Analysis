#!/usr/bin/env python

"""
 mtr_QC.py
 
  
 Built using Anaconda packaged Python:
 

"""

# Standard library.
import datetime, os, sys


# Scientific stack.
import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.dates import MonthLocator, DayLocator, DateFormatter

"""--------------------------------netcdf Routines---------------------------------------"""

def get_global_atts(nchandle):

    g_atts = {}
    att_names = nchandle.ncattrs()
    
    for name in att_names:
        g_atts[name] = nchandle.getncattr(name)
        
    return g_atts

def get_vars(nchandle):
    return nchandle.variables

def ncreadfile_dic(nchandle, params):
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
                data[v] = nchandle.variables[v][:]

        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)
    
"""--------------------------------time Routines---------------------------------------"""

def date2pydate(file_time, file_time2=None, file_flag='EPIC'):


    if file_flag == 'EPIC':
        ref_time_py = datetime.datetime.toordinal(datetime.datetime(1968, 5, 23))
        ref_time_epic = 2440000
    
        offset = ref_time_epic - ref_time_py
    
       
        try: #if input is an array
            python_time = [None] * len(file_time)

            for i, val in enumerate(file_time):
                pyday = file_time[i] - offset 
                pyfrac = file_time2[i] / (1000. * 60. * 60.* 24.) #milliseconds in a day
        
                python_time[i] = (pyday + pyfrac)

        except:
    
            pyday = file_time - offset 
            pyfrac = file_time2 / (1000. * 60. * 60.* 24.) #milliseconds in a day
        
            python_time = (pyday + pyfrac)
        
    else:
        print "time flag not recognized"
        sys.exit()
        
    return np.array(python_time)    
    
"""--------------"""

filein='/Users/bell/Programs/Python/MooringDataProcessing/mtr/pre_corr/15ckpitae_mt4044_0034m.unqcd.orig.nc'
nchandle = Dataset(filein,'a')
global_atts = get_global_atts(nchandle)
vars_dic = get_vars(nchandle)
data = ncreadfile_dic(nchandle, vars_dic.keys())
nchandle.close()

filein='/Users/bell/ecoraid/2015/Moorings/15ckpitae/initial_archive/15ckpitae_s56_0033m.unqcd.nc'
nchandle = Dataset(filein,'a')
global_atts = get_global_atts(nchandle)
vars_dic = get_vars(nchandle)
data2 = ncreadfile_dic(nchandle, vars_dic.keys())
nchandle.close()

filein='/Users/bell/Programs/Python/MooringDataProcessing/mtr/15ckpitae_mt4044_0034m.unqcd.interpolated.trimmed_missing.nc'
nchandle = Dataset(filein,'a')
global_atts = get_global_atts(nchandle)
vars_dic = get_vars(nchandle)
data3 = ncreadfile_dic(nchandle, vars_dic.keys())
nchandle.close()

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
p1 = ax1.plot(date2pydate(data['time'],data['time2']),data['T_20'][:,0,0,0],'b')
date2pydate(data['time'],data['time2'])
ax1.set_ylim([-4,10])
ax1.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
ax1.plot(date2pydate(data2['time'],data2['time2']),data2['T_20'][:,0,0,0],'r')
p1 = ax1.plot(date2pydate(data['time'],data['time2']),data['T_20'][:,0,0,0],'k.')
p1 = ax1.plot(date2pydate(data3['time'],data3['time2']),data3['T_20'][:,0,0,0],'g')