#!/usr/bin/env

"""
 SS_GOA_correlate.py
 
 Correlate two stations (globec and SS3)
 
 Files are created by GOA_Winds_NARR_model_prep.py / SS_Winds_NARR_model_prep.py
 
 -Filtered NARR winds with a triangular filter (1/4, 1/2, 1/4) and output every 3hrs
 -Provided U, V 

 -Saved in EPIC NetCDF standard
"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np

# User Stack
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 04, 29)
__modified__ = datetime.datetime(2014, 04, 29)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR','globec', 'SS3','', 'U,V','Winds', 'Gulf of Alaska'

"""------------------------General   Modules-------------------------------------------"""

def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)
    
def date2pydate(file_time, file_time2=None, file_flag='EPIC'):
    """ Ingest EPIC date or NCEP Date and provide python serial date"""

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
    elif file_flag == 'NCEP':
        """ Hours since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
    elif file_flag == 'Index':
        """ yyyy mm dd"""
        python_time=datetime.datetime.strptime(file_time,'%Y %m %d').toordinal()
    else:
        print "time flag not recognized"
        sys.exit()
        
    return python_time


"""------------------------- MATH   Modules -------------------------------------------"""
def hourly_2_ave(ltbound,utbound, time, data, time_base=6.):
    """ bin average times into specified bins """
    interval = time_base / 24.
    tarray = np.arange(ltbound, utbound,interval)
    dmean = np.zeros_like(tarray) * np.nan
    dstd = np.zeros_like(tarray) * np.nan

    for i, val in enumerate(tarray):
        ind = (time >= val) & (time < val+interval)
        dmean[i] = data[ind].mean()
        dstd[i] = data[ind].std()
    
    return  { 'dtime':tarray, 'dmean':dmean ,'dstd':dstd,} 
    
"""------------------------- Main Routines  -------------------------------------------"""

### NARR wind files (preprocessed) for specific locations - winds have a triangle filter on them
NARR = '/Users/bell/Programs/Python/FOCI_Analysis/'

station_name = ['Globec3','GorePt']
sta_lat = [59.273701,58.9666666666666667]
sta_long = [148.9653,150.9333333333333333]
station_name = [ 'SS2','SS3']
sta_lat = [59.0, 57.6,]
sta_long = [153., 155.]

#loop over all requested data    
years = range(1984, 2014,1)

NARR_globec_time = []
NARR_globec_uwnd = []
NARR_globec_vwnd = []
for iyear in years:
    globec3_data, NARRkeys = from_netcdf(NARR+'GOA_Winds/data/NARR_globec_'+str(iyear)+'.nc')
    NARR_globec_time = NARR_globec_time + date2pydate(globec3_data['time'], globec3_data['time2'])
    NARR_globec_uwnd = np.append(NARR_globec_uwnd, globec3_data['WU_422'][:,0,0,0])
    NARR_globec_vwnd = np.append(NARR_globec_vwnd, globec3_data['WV_423'][:,0,0,0])    

NARR_globec_time = np.array(NARR_globec_time)

### daily averages
time_bin = 24.

NARRDaily_globec_uwnd = hourly_2_ave(NARR_globec_time.min(),NARR_globec_time.max(), NARR_globec_time, NARR_globec_uwnd, time_base=time_bin)    
NARRDaily_globec_vwnd = hourly_2_ave(NARR_globec_time.min(),NARR_globec_time.max(), NARR_globec_time, NARR_globec_vwnd, time_base=time_bin)    

NARR_globec_wndmag = np.sqrt((NARRDaily_globec_uwnd['dmean']**2)+(NARRDaily_globec_vwnd['dmean']**2))

"----"

NARR_ss3_time = []
NARR_ss3_uwnd = []
NARR_ss3_vwnd = []
for iyear in years:
    ss33_data, NARRkeys = from_netcdf(NARR+'SS_Winds/data/NARR_SS3_'+str(iyear)+'.nc')
    NARR_ss3_time = NARR_ss3_time + date2pydate(ss33_data['time'], ss33_data['time2'])
    NARR_ss3_uwnd = np.append(NARR_ss3_uwnd, ss33_data['WU_422'][:,0,0,0])
    NARR_ss3_vwnd = np.append(NARR_ss3_vwnd, ss33_data['WV_423'][:,0,0,0])    

NARR_ss3_time = np.array(NARR_ss3_time)

### daily averages
time_bin = 24.

NARRDaily_ss3_uwnd = hourly_2_ave(NARR_ss3_time.min(),NARR_ss3_time.max(), NARR_ss3_time, NARR_ss3_uwnd, time_base=time_bin)    
NARRDaily_ss3_vwnd = hourly_2_ave(NARR_ss3_time.min(),NARR_ss3_time.max(), NARR_ss3_time, NARR_ss3_vwnd, time_base=time_bin)    

NARR_ss3_wndmag = np.sqrt((NARRDaily_ss3_uwnd['dmean']**2)+(NARRDaily_ss3_vwnd['dmean']**2))

for drange in range(1985,2014,1):
    start_ind = datetime.datetime.strptime(str(drange) + ' 5 01','%Y %m %d').toordinal()
    end_ind = datetime.datetime.strptime(str(drange) + ' 9 15','%Y %m %d').toordinal()
    NARR_ss3_ind = (NARRDaily_ss3_uwnd['dtime'] >= start_ind) & (NARRDaily_ss3_uwnd['dtime'] <= end_ind)
    NARR_globec_ind = (NARRDaily_globec_uwnd['dtime'] >= start_ind) & (NARRDaily_globec_uwnd['dtime'] <= end_ind)
    corr_ss3_globec = np.corrcoef(NARR_globec_wndmag[NARR_globec_ind],NARR_ss3_wndmag[NARR_ss3_ind])
    print "The correlations for May 01 - Sept 15 %s between SS3 and Globec3 are: %s" % (drange, corr_ss3_globec[0][1])
    
    
    