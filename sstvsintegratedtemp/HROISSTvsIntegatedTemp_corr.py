#!/usr/bin/env

"""
 HROISSTvsIntegatedTemp_corr.py
 
 Using NCEP HighRes OI SST and integrated temp - calculate coorelations between the two
 	for each month since 1995

 OISST is Daily and integrated temp is 6hrly so average the integrated temp data
 
"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date

import collections
from scipy import stats

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 06, 06)
__modified__ = datetime.datetime(2016, 06, 06)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NCEP','SST', 'Integrated Temp','M2'


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

def repl_var(nchandle, var_name, val=1e35):
    if len(val) == 1:
        nchandle.variables[var_name][:] = np.ones_like(nchandle.variables[var_name][:]) * val
    else:
        nchandle.variables[var_name][:] = val
    return

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
"""----------------------------- Read from Excel files --------------------------------"""

def readXlsx(fileName,**args):
 import zipfile
 from xml.etree.ElementTree import iterparse
 if "sheet" in args:
    sheet=args["sheet"]
 else:
    sheet=1
 if "header" in args:
    isHeader=args["header"]
 else:
    isHeader=False

 rows = []
 row = {}
 header = {}
 z=zipfile.ZipFile(fileName)
 # Get shared strings
 strings = [el.text for e, el in iterparse(z.open('xl/sharedStrings.xml')) if el.tag.endswith('}t')]
 value = ''

 # Open specified worksheet
 for e, el in iterparse(z.open('xl/worksheets/sheet%d.xml'%(sheet))):
    # get value or index to shared strings
    if el.tag.endswith('}v'): # <v>84</v>
        value = el.text
    if el.tag.endswith('}c'): # <c r="A3" t="s"><v>84</v></c>
        # If value is a shared string, use value as an index
        if el.attrib.get('t') == 's':
            value = strings[int(value)]
        # split the row/col information so that the row leter(s) can be separate
        letter = el.attrib['r'] # AZ22
        while letter[-1].isdigit():
            letter = letter[:-1]
        # if it is the first row, then create a header hash for the names
        # that COULD be used
        if rows ==[]:
            header[letter]=value
        else:
            if value != '': 
                # if there is a header row, use the first row's names as the row hash index
                if isHeader == True and letter in header:
                    row[header[letter]] = value
                else:
                    row[letter] = value

        value = ''
    if el.tag.endswith('}row'):
        rows.append(row)
        row = {}
 z.close()
 return rows

##
# Convert an Excel number (presumed to represent a date, a datetime or a time) into
# a Python datetime.datetime
# @param xldate The Excel number
# @param datemode 0: 1900-based, 1: 1904-based.
# <br>WARNING: when using this function to
# interpret the contents of a workbook, you should pass in the Book.datemode
# attribute of that workbook. Whether
# the workbook has ever been anywhere near a Macintosh is irrelevant.
# @return a datetime.datetime object, to the nearest_second.
# <br>Special case: if 0.0 <= xldate < 1.0, it is assumed to represent a time;
# a datetime.time object will be returned.
# <br>Note: 1904-01-01 is not regarded as a valid date in the datemode 1 system; its "serial number"
# is zero.
# @throws XLDateNegative xldate < 0.00
# @throws XLDateAmbiguous The 1900 leap-year problem (datemode == 0 and 1.0 <= xldate < 61.0)
# @throws XLDateTooLarge Gregorian year 10000 or later
# @throws XLDateBadDatemode datemode arg is neither 0 nor 1
# @throws XLDateError Covers the 4 specific errors

def xldate_as_datetime(xldate, datemode):
    if datemode not in (0, 1):
        sys.exit()
        #raise XLDateBadDatemode(datemode)
    if xldate == 0.00:
        return datetime.time(0, 0, 0)
    if xldate < 0.00:
        sys.exit()
        #raise XLDateNegative(xldate)
    xldays = int(xldate)
    frac = xldate - xldays
    seconds = int(round(frac * 86400.0))
    assert 0 <= seconds <= 86400
    if seconds == 86400:
        seconds = 0
        xldays += 1
    if xldays == 0:
        # second = seconds % 60; minutes = seconds // 60
        minutes, second = divmod(seconds, 60)
        # minute = minutes % 60; hour    = minutes // 60
        hour, minute = divmod(minutes, 60)
        return datetime.time(hour, minute, second)
    return (
        datetime.datetime.fromordinal(xldays + 693594 + 1462 * datemode)
        + datetime.timedelta(seconds=seconds)
        )
    
    

"""--------------------------------main Routines---------------------------------------"""

'''
###Integrated Temperature Prep
int_temp_file = '/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/krafla/ht_anomaly_M2/1995_2015_M2_htcontent_integrated.nc'
nchandle = Dataset(int_temp_file,'a')
global_atts = get_global_atts(nchandle)
vars_dic = get_vars(nchandle)
data1 = ncreadfile_dic(nchandle,vars_dic.keys())
time1 = date2pydate(data1['time'],data1['time2'])
nchandle.close()

plot_var = 'V00_1900'
alltemp = data1[plot_var][:,0,0,0]/74.

datetime1 = num2date(time1,'days since 0001-01-02')

#data is even spaced but starts on 6am of 3/17/95
subtime=[]
for now in datetime1:
    if now.hour == 0:
        subtime = subtime + [now]

t=alltemp[3:].reshape(7594,4)
alltemp_daily=t.mean(axis=1)
alltemp_daily[np.where(alltemp_daily>1e30)] = 1e35
    
for i,v in enumerate(subtime):
    print v, alltemp_daily[i]
'''
''' 
monthlymeans = np.array([])

###SST Files
path_sst = '/Volumes/WDC_internal/Users/bell/Data_Local/sst/M2_NOAA_OI_SST_V2/NOAA_OI_SST_V2_M2_'
ssttime, sst = [], []
for year in range(1995, 2016,1):
    nchandle = Dataset(path_sst + str(year) + '.nc','a')
    global_atts = get_global_atts(nchandle)
    vars_dic = get_vars(nchandle)
    data = ncreadfile_dic(nchandle,vars_dic.keys())
    ssttime = np.hstack((ssttime, date2pydate(data['time'],data['time2'])))
    sst = np.hstack((sst, data['T_25'][:,0,0,0]))
    nchandle.close()

ssttime = num2date(ssttime,'days since 0001-01-02')
'''

### read in data file
# from each line, look for expected values to be archived filling in missing
# values with 1e35.  Build nc files as each cast cycles
excelfiles = '/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/stabeno/june/sstvsinttemp/inttemp_oisst_correlationtable.xlsx'
print "Reading file {0}".format(excelfiles)
W = readXlsx(excelfiles, sheet=2, header=True)
data = collections.OrderedDict()
for index, row in enumerate(W):
    if index<=1-1:
        print row
    else:
        try:
            data[index] = row
        except:
            print "End of data stream"    

print "yyyymm, samplnum, sstmean, integrated temp mean, rcoef, rcoefsqrd"
for year in range(1995,2016,1):
    for month in range(1,13,1):
        ind = []
        x,y = [], []
        for k in data.keys():
            tempdate = xldate_as_datetime(float(data[k]['date']),0)

            if (tempdate.year == year) and (tempdate.month == month):
                ind = ind + [k]
                if not float(data[k]['int_temp']) > 1e30: #skip missing data
                    x = x + [float(data[k]['oi_sst'])]
                    y = y + [float(data[k]['int_temp'])]
        try:        
            lstats = stats.linregress(x,y)
            print "{0}-{1}, {2}, {3}, {4}, {5}, {6}".format(year, month, len(y), np.array(x).mean(), np.array(y).mean(), lstats[2], lstats[2]**2.)
        except:
            print "{0}-{1}, {2}, {3}, {4}, {5}, {6}".format(year, month, 'na','na','na','na', 'na')
            