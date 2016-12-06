#!/usr/bin/env

"""
Program: nc2r_cluster.py

Purpose: take a netcdf 3dim file and subset it based on time, space

Input/Arguments: Latitude and Longitude, raw data file
 
"""
#System Stack
import datetime

#Science Stack
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np

# User Stack
import utilities.ncutilities as ncutil


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 03, 25)
__modified__ = datetime.datetime(2014, 03, 25)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'R analysis', 'SST','Chukchi'

"""--------------------------------netcdf Routines---------------------------------------"""


def get_global_atts(nchandle):

    g_atts = {}
    att_names = nchandle.ncattrs()
    
    for name in att_names:
        g_atts[name] = nchandle.getncattr(name)
        
    return g_atts

def get_vars(nchandle):
    return nchandle.variables

def get_var_atts(nchandle, var_name):
    return nchandle.variables[var_name]

def ncreadfile_dic(nchandle, params):
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
                data[v] = nchandle.variables[v][:]

        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)

def ncreadfile_dic_oneday(nchandle, params,doy,coords):
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
            if v in coords:
                data[v] = nchandle.variables[v][:]
            else:
                data[v] = nchandle.variables[v][doy]

        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)
"""--------------------------------time Routines---------------------------------------"""

def date2pydate(file_time, file_time2=None, time_since_str=None, file_flag='EPIC'):


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
        

    elif file_flag == 'NARR':
        """ Hours since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
    elif file_flag == 'NCEP_days':
        """ days since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time + base_date
    elif file_flag == 'NCEP':
        """ Hours since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
    elif file_flag == 'netCDF4':
        """ Use time_since_str"""
        python_time = num2date(file_time,'time_since_str','standard')
    else:
        print "time flag not recognized"
        sys.exit()
        
    return np.array(python_time)

"""------------------------------------- MAPS -----------------------------------------"""



def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx

"""------------------------------- MAIN--------------------------------------------"""

ncfile_all = ['/Users/bell/in_and_outbox/2016/bond/feb/ssh_a.nc', ]

for ncfile in ncfile_all:
    print "Working with {0}".format(ncfile)
    nchandle = Dataset(ncfile,'r')
    params = get_vars(nchandle)
    data = ncreadfile_dic(nchandle, params)
    nchandle.close()

lon_bounds = [180,203] #degrees_east
lat_bounds = [66,72] #degrees North
time_bounds = [6, 7, 8, 9, 10] #june-nov

lon_ind = np.where((data['LON171_220'] >= lon_bounds[0]) & (data['LON171_220'] <= lon_bounds[1]) )[0]
lat_ind = np.where((data['LAT151_165'] >= lat_bounds[0]) & (data['LAT151_165'] <= lat_bounds[1]) )[0]

#print header
for lon in data['LON171_220'][lon_ind]:
    for lat in data['LAT151_165'][lat_ind]:
        print "lo{0}la{1}".format(lon,lat)
    
#find valid times, print each time 
for i,v in enumerate(data['TMON']):
    temptime = num2date(v,'DAYS since 1901-01-15 00:00:00','standard')
    if temptime.month in time_bounds:
        temp_data = ''
        for lo in lon_ind:
            for la in lat_ind:
                print temp_data
                temp_data = temp_data + ', ' + str(data['SSH_A'].data[i][la][lo])
        print 't' + str(temptime.year) + '_' + str(temptime.month) + temp_data
        