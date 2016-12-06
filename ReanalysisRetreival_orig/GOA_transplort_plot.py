#!/usr/bin/env python

"""
 GOA_transplort_plot.py
 
 Plot transport data from GOA for 2015 ngoa paper.
 
 Reads EPIC netcdf of transport as a function of time

"""

# Standard library.
import datetime, sys

# System Stack
import argparse

# Scientific stack.
import numpy as np
from netCDF4 import Dataset

# Visual Stack
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2015, 06, 17)
__modified__ = datetime.datetime(2015, 06, 17)
__version__  = "0.1.0"
__status__   = "Development"

        
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
    
"""------------------------------- MAIN ----------------------------------------"""

parser = argparse.ArgumentParser(description='EPIC NetCDF Timeseries Plot')
parser.add_argument('DataPath', metavar='DataPath', type=str, help='full path to file')
parser.add_argument("-monthly_average",'--monthly_average', action="store_true", help='calculate monthly average')
parser.add_argument("-annual_signal",'--annual_signal', action="store_true", help='calculate total annual signal')

args = parser.parse_args()

###nc readin
nchandle = Dataset(args.DataPath,'a')
global_atts = get_global_atts(nchandle)
vars_dic = get_vars(nchandle)
data1 = ncreadfile_dic(nchandle, vars_dic.keys())

nchandle.close()

### convert epic time to python serialdate
time1 = date2pydate(data1['time'],data1['time2'])
xdata = data1['TR_388'][:,0,0,0]
label = 'TR_388'

scale = 10**6

if args.monthly_average or args.annual_signal:

    proc_data = {}
    ave_data = {}
    ave_data_stats = {}
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    for i,v in enumerate(time1):
        proc_data[datetime.datetime.fromordinal(int(v)).strftime('%Y-%m-%d')] = data1['TR_388'][i,0,0,0]
    
    for years in np.arange(1984,2007):
        for month in months:
            sum = 0
            count = 0
            for k in proc_data.keys():
                if str(years)+'-'+month in k:
                    if proc_data[k] < 1e30:
                        count += 1 
                        sum = sum + proc_data[k]
            print str(years)+'-'+month+'-15'
            if count != 0:
                ave_data[str(years)+'-'+month+'-15'] = sum / scale / count
                ave_data_stats[str(years)+'-'+month+'-15'] = count
            else:
                ave_data[str(years)+'-'+month+'-15'] = 1e35
                ave_data_stats[str(years)+'-'+month+'-15'] = 0
    
    
    time1 = []
    xdata = []
    
    ### build x and y data
    # for x data we want all the months of the entire record to show up at the same 'x-point' so every feb for all years is collocated
    for k in sorted(ave_data.keys()):
        time1 = time1 + [datetime.datetime.strptime("2000-"+"-".join(k.split('-')[1:]),'%Y-%m-%d').toordinal()]
        xdata = xdata + [ave_data[k]]
    time1 = np.array(time1)
    xdata = np.array(xdata)
        
    if args.annual_signal:
        timet = []
        xdatat = []
        for month in months:
            time_base = datetime.datetime.strptime("2000-"+month+"-15",'%Y-%m-%d').toordinal()
            ave_ind = np.where(time1==time_base) #arbitrary 2000 year chosen above for plotting
            dave_ind = np.where(xdata[ave_ind] != 1e35)
            timet = timet + [time_base]
            xdatat = xdatat + [np.mean(xdata[ave_ind][dave_ind])]
        time1 = np.array(timet)
        xdata = np.array(xdatat)

### plot variable against time
fig = plt.figure(2)

ax2 = plt.subplot2grid((3, 1), (1, 0), colspan=1, rowspan=3)
p2 = ax2.plot(time1, xdata,'r.', markersize=12)
#ax2.set_ylim([xdata[xdata != 1e35].min()-0.5,xdata[xdata != 1e35].max()+0.5])
ax2.set_ylim([-2,3])
ax2.set_xlim([datetime.datetime.strptime("2000-1-1",'%Y-%m-%d').toordinal(),datetime.datetime.strptime("2001-1-1",'%Y-%m-%d').toordinal()])
plt.ylabel(label)
ax2.xaxis.set_major_locator(MonthLocator(interval=1))
ax2.xaxis.set_minor_locator(MonthLocator())
ax2.xaxis.set_major_formatter(DateFormatter('%b'))

fig.autofmt_xdate()

fullpath = '/Users/bell/temp/'
user_append = 'minus_downstream_monthlyaverage_ma'
plt.savefig(fullpath + user_append + '_ngoa_fig8.svg', bbox_inches='tight', dpi = (300))

plt.close()
