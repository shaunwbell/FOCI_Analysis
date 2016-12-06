#!/usr/bin/env

"""
 epic_make_monthly.py
 
 From a list of netcdf files calculate for chosen parameter
    30.42 day average and output
    
    year    month   var_ave number_of_samples
    2001    1   30.1    15
    2001    2   30.0    14
    2001    3   1e35    0    
    

    print output to screen to be captured

Example:
--------
python epic_make_monthly.py cb1_uv_brg_brg_f35.point U_320 > cb1_mean_U320.txt

"""
#System Stack
import datetime
import argparse
from netCDF4 import Dataset

#Science Stack
import numpy as np


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 04, 29)
__modified__ = datetime.datetime(2014, 04, 29)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'monthly','averages'


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
    
"""------------------------------- MAIN--------------------------------------------"""

parser = argparse.ArgumentParser(description='make 30.42 day monthly files from timeseries data')
parser.add_argument('inputpath', metavar='inputpath', type=str, help='full path pointer file with paths to data')
#parser.add_argument('output', metavar='output', type=str, help='optional output path')
parser.add_argument('epic_key', metavar='epic_key', type=str, help='epic key code for variable to be processed')
args = parser.parse_args()

#ptr file has one mooring per line
if '.txt' in args.inputpath:  
    with open(args.inputpath,'r') as fid:
        ifile = fid.read()
    fid.close()
    nc_files = ifile.split('\n')
elif '.point' in args.inputpath:  
    with open(args.inputpath,'r') as fid:
        ifile = fid.read()
    fid.close()
    nc_files = ifile.strip().split('\n')
else:
    nc_files = [args.inputpath,]
        
### READ .nc files -- assumes the order of files in ptr file is increasing time
## Only retain designated variable
time = []
var = []
for ncfile in nc_files:
    print "Reading file {0}".format(ncfile)
    nchandle = Dataset(ncfile,'r')
    global_atts = get_global_atts(nchandle)
    vars_dic = get_vars(nchandle)
    data = ncreadfile_dic(nchandle, vars_dic.keys())
    nchandle.close()
    
    time = np.hstack((time,date2pydate(data['time'],data['time2'])))
    var = np.hstack((var,data[args.epic_key][:,0,0,0]))
    

### After all files are read in, find smallest date
year_min = datetime.datetime.fromordinal(int(np.min(time))).year
year_max = np.ceil(np.max(time))
year_ord = datetime.datetime.toordinal(datetime.datetime(year_min,1,1))

month = 1
window = 30.42
for mind in np.arange(year_ord,year_max,window):
    data_ind = np.where((time >= mind) & (time < mind + window))
    missing_ind = np.where(var[data_ind] < 1e30)    
    window_average = np.average(var[data_ind][missing_ind])
    if np.isnan(window_average):
        window_average = 1e35
    
    print "{0} {1} {2} {3}".format(datetime.datetime.fromordinal(int(np.min(mind))).year, month,
        window_average, len(var[data_ind][missing_ind]) )
        
    if month == 12:
        month = 1
    else:
        month +=1
