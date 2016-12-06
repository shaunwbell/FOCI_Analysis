#!/usr/bin/env

"""
 NCEP_NARR_comp.py
 
 NCEP vs NARR side by side comparisons of select fields
 
 Compare NARR Winds with NCEP V2 (with Mooring Winds)

 Using Anaconda packaged Python 
"""
#System Stack
import datetime

#Science Stack
import numpy as np

# User Stack
import general_utilities.date2doy as date2doy
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.dates import MonthLocator, DateFormatter

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"

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
        base_date=datetime.datetime.strptime(file_time,'%Y-%m-%d').toordinal()
        python_time = file_time * 24. + base_date
    elif file_flag == 'NCEP':
        """ Hours since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
    else:
        print "time flag not recognized"
        sys.exit()
        
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
    
def cart2wind(cart_angle):
    """ 0deg is North, rotate clockwise"""
    cart_angle_out = np.zeros_like(cart_angle)
    cart_angle = 90. - cart_angle  #rotate so N is 0deg
    cart_angle =cart_angle % 360.
    
    return cart_angle
    
def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)

    
"""---------------------------- Main Routine-------------------------------------------"""

"""------Ingest Data--------"""

#doy = date2doy.date2doy('2003-06-03')

NCEPV2 = '/Users/bell/Data_Local/Reanalysis_Files/NCEPV2/daily_mean/'

NCEPV2_uwind, NCEPV2_uparams = from_netcdf(NCEPV2 + 'uwnd.10m.gauss.2003.nc')
NCEPV2_vwind, NCEPV2_vparams = from_netcdf(NCEPV2 + 'vwnd.10m.gauss.2003.nc')
NCEPTime = date2pydate(NCEPV2_uwind['time'], file_flag='NCEP')

NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'

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
MooringDaily_uwnd = hourly2daily(NARRTime.min(),NARRTime.max(), MooringTime, MooringMetData['WU_422'])
MooringDaily_vwnd = hourly2daily(NARRTime.min(),NARRTime.max(), MooringTime, MooringMetData['WV_423'])
sta_lat = MooringMetData['latitude'][0]
sta_long = MooringMetData['longitude'][0]

"""---------------------------- Data Manipulation Routines-----------------------------"""


#NCEP V2 - vectors given to define grid
#shift grid from 0->360 to -360->0
NCEP_uwnds,lons_ncep = shiftgrid(0.,NCEPV2_uwind['uwnd'],NCEPV2_uwind['lon'],start=False) 
NCEP_vwnds,lons_ncep = shiftgrid(0.,NCEPV2_vwind['vwnd'],NCEPV2_vwind['lon'],start=False) 

if isinstance(NARR_uwind['uwnd'], np.ma.core.MaskedArray): #masked array handling
    NARR_uwind['uwnd'] = NARR_uwind['uwnd'].data
    NARR_vwind['vwnd'] = NARR_vwind['vwnd'].data
    
### Plot

for ind in (range(0,366,1)):

    # use for color coding of winds
    CmagNARR = np.sqrt(NARR_uwind['uwnd'][ind,:,:]**2. + NARR_vwind['vwnd'][ind,:,:]**2.)
    CmagNCEP = np.sqrt(NCEP_uwnds[ind,:,:][0]**2. + NCEP_vwnds[ind,:,:][0]**2.)

    fig = plt.figure()
    
    ax = plt.subplot(121)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=55, \
        urcrnrlat=65,llcrnrlon=-155,urcrnrlon=-140, lat_ts=45)


    #NARR - array given to define grid
    x_narr, y_narr = m(NARR_uwind['lon'], NARR_uwind['lat'])

    # Mooring Data
    x_moor, y_moor = m(-1. * sta_long,sta_lat)
    Q1 = m.quiver(x_moor,y_moor,MooringDaily_uwnd['daily_mean'][ind],MooringDaily_vwnd['daily_mean'][ind],scale=100, color='b')

    Q = m.quiver(x_narr,y_narr,NARR_uwind['uwnd'][ind,:,:],NARR_vwind['vwnd'][ind,:,:], CmagNARR, cmap='Reds', scale=100)
    qk = plt.quiverkey(Q, 0.05, 0.05, 5, '5 m/s', labelpos='S')

    m.scatter(x_moor,y_moor,10,marker='o',color='g')
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(55,65,2.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-155,-140,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    ax = plt.subplot(122)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=55, \
        urcrnrlat=65,llcrnrlon=-155,urcrnrlon=-140, lat_ts=45)

    lon_ncep, lat_ncep = np.meshgrid(lons_ncep, NCEPV2_uwind['lat'])
    x_ncep, y_ncep = m(lon_ncep, lat_ncep)


    # Mooring Data
    x_moor, y_moor = m(-1. * sta_long,sta_lat)
    Q1 = m.quiver(x_moor,y_moor,MooringDaily_uwnd['daily_mean'][ind],MooringDaily_vwnd['daily_mean'][ind],scale=100, color='b')
    
    Q = m.quiver(x_ncep,y_ncep,NCEP_uwnds[ind,:,:][0],NCEP_vwnds[ind,:,:][0], CmagNCEP, cmap='Reds', scale=100)
    
    #
    qk = plt.quiverkey(Q, 0.05, 0.05, 5, '5 m/s', labelpos='S')

    m.scatter(x_moor,y_moor,10,marker='o',color='g')
    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(55,65,2.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-155,-140,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    ind_off = ind + 1 #offset array with DOY
    if (ind_off < 10):
        str_ind = '00' + str(ind_off)
    elif (ind_off >= 10 and ind_off < 100):
        str_ind = '0' + str(ind_off)
    else:
        str_ind = str(ind_off)
    
    fig.suptitle('NARR (left) vs NCEP V2 (right) DOY:'+str_ind, fontsize=12)
    
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*1.25, DefaultSize[1]) )
    plt.savefig('images_2003/Globec_region' + str_ind + '.png', bbox_inches='tight', dpi = (100))
    plt.close(fig)

   
    print "Finishing Figure: " + str_ind
