#!/usr/bin/env

"""
 GOA_Winds_NARR_3hr.py
 
 Compare NARR Winds with NCEP V2 (with Mooring Winds)

 Using Anaconda packaged Python 
"""
#System Stack
import datetime

#Science Stack
import numpy as np

# User Stack
import general_utilities.date2doy as date2doy
import general_utilities.haversine as sphered
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.dates import MonthLocator, DateFormatter, DayLocator

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR','GLOBEC3','3hr comparison', 'Winds', 'Gulf of Alaska'

"""------------------------General   Modules-------------------------------------------"""

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
            pytime = []
    
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
    else:
        print "time flag not recognized"
        
    return np.array(python_time)
        
def hourly2daily(ltbound,utbound, time, data):
    tarray = np.arange(ltbound, utbound+1,1.)
    dmean = np.zeros_like(tarray)
    dstd = np.zeros_like(tarray)
    
    for i, val in enumerate(tarray):
        ind = np.where(np.floor(time) == val )
        dmean[i] = data[ind].mean()
        dstd[i] = data[ind].std()
    

    return ( {'daily_mean':dmean ,'daily_std':dstd, 'daily_time':tarray} )

def hourly_2_3hrly(ltbound,utbound, time, data):
    interval = 3 / 24.
    tarray = np.arange(ltbound, utbound,interval)
    dmean = np.zeros_like(tarray)
    dstd = np.zeros_like(tarray)
    
    for i, val in enumerate(tarray):
        ind = (time >= val) & (time < val+interval)
        dmean[i] = data[ind].mean()
        dstd[i] = data[ind].std()
    

    return ( {'mean':dmean ,'std':dstd, 'time':tarray} )
    
        
def cart2wind(cart_angle):
    """ 0deg is North, rotate clockwise"""

    cart_angle = 90. - cart_angle  #rotate so N is 0deg
    cart_angle =cart_angle % 360.
    
    return cart_angle
    
def rotate_coord(angle_rot, mag, dir):
    """ converts math coords to along/cross shelf.
    + onshore  / along coast with land to right (right handed)
    - offshore / along coast with land to left
    
    Todo: convert met standard for winds (left handed coordinate system
    """
    
    dir =  dir - angle_rot
    
    along = mag * np.sin(np.deg2rad(dir))
    cross = mag * np.cos(np.deg2rad(dir))
    
    return (along, cross)
    
def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)

"""---------------------------- Plotting Modules --------------------------------------"""

def quiver_timeseries(time,ucomp,vcomp,magnitude,data_source):
    
    t_ind = ~(~np.isnan(magnitude) & (magnitude < 100))
    ucomp[t_ind] = 0.
    vcomp[t_ind] = 0.
    magnitude[t_ind] = 0.
        
    fig1, (ax1, ax2) = plt.subplots(2,1)
    # Plot quiver
    ax1.set_ylim(-magnitude.max(), magnitude.max())
    fill1 = ax1.fill_between(time, magnitude, 0, color='k', alpha=0.1)

    # Fake 'box' to be able to insert a legend for 'Magnitude'
    p = ax1.add_patch(plt.Rectangle((1,1),1,1,fc='k',alpha=0.1))
    leg1 = ax1.legend([p], ["Wind magnitude [m/s]"],loc='lower right')
    leg1._drawFrame=False

    # 1D Quiver plot
    q = ax1.quiver(time,0,ucomp,vcomp,color='r',units='y',scale_units='y',
                   scale = 1,headlength=1,headaxislength=1,width=0.04,alpha=.95)
    qk = plt.quiverkey(q,0.2, 0.05, 5,r'$5 \frac{m}{s}$',labelpos='W',
                   fontproperties={'weight': 'bold'})

    # Plot u and v components
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_xlim(time.min(),time.max()+0.5)
    ax1.set_ylabel("Velocity (m/s)")
    ax2.plot(time, vcomp, 'b-')
    ax2.plot(time, ucomp, 'g-')
    ax2.set_xlim(time.min(),time.max()+0.5)
    ax2.set_xlabel("Date (UTC)")
    ax2.set_ylabel("Velocity (m/s)")
    ax2.xaxis.set_major_locator(MonthLocator())
    ax2.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    ax2.xaxis.set_minor_locator(DayLocator())
    fig1.autofmt_xdate()    
    # Set legend location - See: http://matplotlib.org/users/legend_guide.html#legend-location
    leg2 = plt.legend(['v','u'],loc='upper left')
    leg2._drawFrame=False
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('both')

    
    DefaultSize = fig1.get_size_inches()
    fig1.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )

    fig1.suptitle("3hr ave Wind data for:   " + data_source, fontsize=12)
    # Save figure (without 'white' borders)
    plt.savefig('images/Globec_' + data_source + '_timeseries.png', bbox_inches='tight', dpi = (100))
    plt.close(fig1) 
    
    
"""---------------------------- Main Routine-------------------------------------------"""

"""------Ingest Data--------"""

NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/3hourly/'

NARR_uwind, NARR_uparams = from_netcdf(NARR + 'uwnd.10m.2003.nc')
NARR_vwind, NARR_vparams = from_netcdf(NARR + 'vwnd.10m.2003.nc')
NARRTime = date2pydate(NARR_uwind['time'], file_flag='NCEP')

### NARR Data has the following boundary corners:
# 12.2N;133.5W, 54.5N; 152.9W, 57.3N; 49.4W ,14.3N;65.1W
# Lambert Conformal

#lat/lon is ~ 59N, 149W
MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/2003/globec3/03gbm3a_wpak.nc'
MooringMetData, Mooring_params = from_netcdf(MooringFile)

MooringTime = date2pydate(MooringMetData['time'], MooringMetData['time2'], file_flag='EPIC')
MooringDaily_uwnd = hourly_2_3hrly(NARRTime.min(),NARRTime.max(), MooringTime, MooringMetData['WU_422'])
MooringDaily_vwnd = hourly_2_3hrly(NARRTime.min(),NARRTime.max(), MooringTime, MooringMetData['WV_423'])
sta_lat = MooringMetData['latitude'][0]
sta_long = MooringMetData['longitude'][0]

#-----> user set to force mooring location instead of using built in location (useful if you want
# to specify lat/lon for model comparison purposes
#sta_lat = 58.
#sta_long = 148.
"""---------------------------- Data Manipulation Routines-----------------------------"""

#Find NCEP and NARR nearest point to mooring
narrpt = sphered.nearest_point([sta_lat,-1 * sta_long],NARR_uwind['lat'],NARR_uwind['lon'], '2d')

try: #some data gives extra paramter and puts data in a dict structure... others does not
    NARR_u = NARR_uwind['uwnd'].data[:,narrpt[3],narrpt[4]]
    NARR_v = NARR_vwind['vwnd'].data[:,narrpt[3],narrpt[4]]
    
except TypeError: #no .data parameter
    NARR_u = NARR_uwind['uwnd'][:,narrpt[3],narrpt[4]]
    NARR_v = NARR_vwind['vwnd'][:,narrpt[3],narrpt[4]]
    
NARR_wind_mag = np.sqrt(NARR_u**2. + NARR_u**2.)
NARR_wind_dir_math = np.rad2deg(np.arctan2(NARR_v , NARR_u)) 
NARR_wind_dir = cart2wind(NARR_wind_dir_math)

Mooring_wind_mag = np.sqrt(MooringDaily_uwnd['mean']**2. + MooringDaily_vwnd['mean']**2.)
Mooring_wind_dir_math = np.rad2deg(np.arctan2(MooringDaily_vwnd['mean'] , MooringDaily_uwnd['mean'])) 
Mooring_wind_dir = cart2wind(Mooring_wind_dir_math)

# mask when mooring wasn't available
t_ind = ~np.isnan(Mooring_wind_mag) & (Mooring_wind_mag < 100)

### Calculate +-flow and x-flow rotating along coast (~43 degrees bearing near Globec3 )
(NARRalong, NARRcross) = rotate_coord(137., NARR_wind_mag, NARR_wind_dir_math)
(MOORalong, MOORcross) = rotate_coord(137., Mooring_wind_mag, Mooring_wind_dir_math)

"""---------------------------- Plotting Routines--------------------------------------"""

### standard wind / time plots
#   NARR
quiver_timeseries(NARRTime,NARR_u,NARR_v,NARR_wind_mag,'NARR')
quiver_timeseries(MooringDaily_uwnd['time'],MooringDaily_uwnd['mean'],MooringDaily_vwnd['mean'],Mooring_wind_mag,'GLOBEC3')

### Along/Cross Shore comparisons Mooring vs NARR/NCEP
#  for entire year (mark mooring specific times)
fig = plt.figure(6)

ax = plt.subplot(121)
p1 = ax.plot(MOORalong[t_ind], NARRalong[t_ind])
plt.setp(p1,color='r', marker='+', linestyle='')
p3 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(-15,16,5))
ax.set_yticks(np.arange(-15,16,5))
ax.set_xlabel('3hr Globec3 Alongshore Flow (m/s)')
ax.set_ylabel('3hr NARR Alongshore Flow (m/s)')

ax = plt.subplot(122)
p1 = ax.plot(MOORcross[t_ind], NARRcross[t_ind])
plt.setp(p1,color='r', marker='+', linestyle='')
p3 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(-15,16,5))
ax.set_yticks(np.arange(-15,16,5))
ax.set_xlabel('3hr Globec3 AcrossShore Flow (m/s)')
ax.set_ylabel('3hr NARR AcrossShore Flow (m/s)')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_alongacross_comp.png', bbox_inches='tight', dpi = (100))
plt.close()
