#!/usr/bin/env

"""
 SS_Winds_HovmellerIMG.py
 
 Hovmoller diagram of wind speed vs time (as an imgsec plot)
    GorePoint - 58deg 58min N, 150deg 56min W 
    and Globec3 59.273701N, 148.9653W

 Files are created by GOA_Winds_NARR_model_prep.py
 
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
__keywords__ = 'NARR','Shelikof', 'SS3','AO/NAO/PNA', 'U,V','Winds', 'Gulf of Alaska'

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
NARR = '/Users/bell/Programs/Python/FOCI_Analysis/SS_Winds/data/'

station_name = [ 'SS2','SS3']
sta_lat = [59.0, 57.6,]
sta_long = [153., 155.]

#loop over all requested data    
years = range(1984, 2014,1)

NARR_time = []
NARR_uwnd = []
NARR_vwnd = []
for iyear in years:
    ss3_data, NARRkeys = from_netcdf(NARR+'NARR_SS3_'+str(iyear)+'.nc')
    NARR_time = NARR_time + date2pydate(ss3_data['time'], ss3_data['time2'])
    NARR_uwnd = np.append(NARR_uwnd, ss3_data['WU_422'][:,0,0,0])
    NARR_vwnd = np.append(NARR_vwnd, ss3_data['WV_423'][:,0,0,0])    

NARR_time = np.array(NARR_time)

### daily averages
time_bin = 24.

NARRDaily_uwnd = hourly_2_ave(NARR_time.min(),NARR_time.max(), NARR_time, NARR_uwnd, time_base=time_bin)    
NARRDaily_vwnd = hourly_2_ave(NARR_time.min(),NARR_time.max(), NARR_time, NARR_vwnd, time_base=time_bin)    

NARR_wndmag = np.sqrt((NARRDaily_uwnd['dmean']**2)+(NARRDaily_vwnd['dmean']**2))

#convert timeseries into an array with x-axis -> April 1 - Oct 1 and y axis being years
y_imNARR = range(1985,2014,1)
x_imNARR = range(datetime.datetime.strptime('1985 4 01','%Y %m %d').toordinal(), 
                 datetime.datetime.strptime('1985 10 01','%Y %m %d').toordinal(), 1)
for drange in y_imNARR:
    start_ind = datetime.datetime.strptime(str(drange) + ' 4 01','%Y %m %d').toordinal()
    end_ind = datetime.datetime.strptime(str(drange) + ' 10 01','%Y %m %d').toordinal()
    NARR_ind = (NARRDaily_uwnd['dtime'] >= start_ind) & (NARRDaily_uwnd['dtime'] <= end_ind)
    if drange == np.min(y_imNARR):
        NARR_mag_array = NARR_wndmag[NARR_ind]
    else:
        NARR_mag_array = np.vstack((NARR_mag_array,NARR_wndmag[NARR_ind]))

fig = plt.figure()
ax = plt.subplot(111)
plt.imshow(NARR_mag_array,cmap='gray_r', aspect=20, 
    origin='lower',extent=[np.min(x_imNARR),np.max(x_imNARR),np.min(y_imNARR),np.max(y_imNARR)])
ax.xaxis.set_major_formatter(DateFormatter('%b %d'))
ax.xaxis.set_major_locator(MonthLocator(bymonth=range(4,10,2),bymonthday=15))
ax.xaxis.set_minor_locator(MonthLocator(bymonthday=[1,15]))
plt.minorticks_on()
cbar = plt.colorbar()

plt.title('NARR Wind Magnitude at SS3 (exit) \n')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )
plt.savefig('NARR_SS3_windmag.png', bbox_inches='tight', dpi = (100))
plt.close()