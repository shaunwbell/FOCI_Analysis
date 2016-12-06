#!/usr/bin/env

"""
 BS_Winds_NARR_6hr.py
 
 Compare NARR Winds with NCEP V2 (with Mooring Winds) for 6hr intervals.  Uses 3hr NARR and 6hr NCEP

 Using Anaconda packaged Python 
"""
#System Stack
import datetime

#Science Stack
import numpy as np
from scipy import stats
from netCDF4 import Dataset

# User Stack
import general_utilities.haversine as sphered
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.dates import MonthLocator, DateFormatter, DayLocator

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 03, 25)
__modified__ = datetime.datetime(2014, 03, 25)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR', 'NCEP V2','Bering2', 'power law', 'user defined time comparison', 'Winds', 'Bering Sea'

"""------------------------General   Modules-------------------------------------------"""

def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)

def from_netcdf_mf(infiles):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.mf_ncopen(infiles)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)
    
def from_netcdf_1dsplice(infile, height_ind, lat_ind, lon_ind):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them

    print "Parameters available: " 
    print params
    
    ncdata = ncutil.ncreadfile_dic_slice(nchandle, params, height_ind=height_ind, lat_ind=lat_ind, lon_ind=lon_ind)
    ncutil.ncclose(nchandle)
    
    return ncdata
    
def missing2nan(MooringMetData, Mooring_params, fill_value=np.nan):
    """replace known placeholder values with nan or other fill_value"""
    for i,v in enumerate(Mooring_params):
        not_nan_ind = (MooringMetData[v] == 1e35)
        try:
            MooringMetData[v][not_nan_ind] = np.nan
        except ValueError:
            pass #no missing data so not_nan_ind is empty
        
    return MooringMetData

def latlon_grid(infile):
    nchandle = ncutil.ncopen(infile)
    lat_lon = ncutil.get_geocoords(nchandle)
    ncutil.ncclose(nchandle)

    return (lat_lon)

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
    elif file_flag == 'NCEPV1':
        """ Hours since 1-1-1"""
        base_date=datetime.datetime.strptime('0001-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
        python_time = python_time - 1.75 #<-- hack correction accounts for python date discrepancy? and 6hr lead/lag of dataset
    else:
        print "time flag not recognized"
        sys.exit()
        
    return np.array(python_time)

def etopo5_data():
    """ read in etopo5 topography/bathymetry. """
    file = '/Users/bell/Data_Local/MapGrids/etopo5.nc'
    etopodata = Dataset(file)
    
    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()
    
    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0
    
    lons, lats = np.meshgrid(lons, lats)
    
    return(topoin, lats, lons)
            
"""------------------------------- Stats/Math Modules --------------------------------------"""

def lin_fit(x, y):
    """ rely's on scipy"""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    return ( slope, intercept, r_value, p_value, std_err )

def comp_corr( x, y):
    """
    Complex Correlations
    
    Parameters:
    -----------
    x: complex vector 1
    y: complex vector 2
    
    Outputs:
    --------
    
    complex correlation vector between x and y (orientation independent)
    complex correlation angle (ccw rotation of y with respect to x)
    
    Reference:
    ----------
    Kundu, Pijush K., 1976: Ekman Veering Observed near the Ocean Bottom. J. Phys. Oceanogr., 6, 238-242
    """
    
    x = x[0] + 1j* x[1]
    y = y[0] + 1j* y[1]
    
    # From equation 3.3
    corr = np.inner(np.conjugate(x),y) \
    / (np.sqrt(np.inner(np.conjugate(x),x)) * np.sqrt(np.inner(np.conjugate(y),y)))
    
    corr_mag = np.sqrt(corr.real**2 +corr.imag**2)
    
    corr_angle = np.rad2deg(np.arctan2(corr.imag, corr.real))
    """
    # From equation 3.6 and 3.7
    
    # what is the innerproduct of <u1u2 + v1v2> ???
    real_c = (x[0]*y[0] + x[1]*y[1]) / (np.sqrt(x[0]**2. + y[0]**2.) * np.sqrt(x[1]**2. + y[1]**2.))
    imag_c = 1j * (x[0]*y[1] - x[1]*y[0]) / (np.sqrt(x[0]**2. + y[0]**2.) * np.sqrt(x[1]**2. + y[1]**2.))

    corr_angle = np.arctan2((x[0]*y[1] - x[1]*y[0]), (x[0]*y[0] + x[1]*y[1]))
    """
    return (corr_mag, corr_angle)
    
def wind_power_law(comp_orig, height_obs=3., height_interp=10., correction=False):
    """simple power law wind adjustment
    default - 3m observations, 10m interpolated height"""
    
    if correction:
        wind_cor = comp_orig * (height_interp / height_obs)**(0.143)
    else:
        wind_cor = comp_orig
        
    return wind_cor
    
def hourly_2_ave(ltbound,utbound, time, data, time_base=6.):
    interval = time_base / 24.
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
    

"""---------------------------- Plotting Modules --------------------------------------"""

def quiver_timeseries(time,ucomp,vcomp,magnitude,data_source,station_name):
    
    t_ind = ~(~np.isnan(magnitude) & (magnitude < 100))
    ucomp[t_ind] = 0.
    vcomp[t_ind] = 0.
    magnitude[t_ind] = 0.
        
    fig1, (ax1, ax2) = plt.subplots(2,1)
    # Plot quiver
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
    # Plot u and v components
    #ax1.set_ylim(-magnitude.max(), magnitude.max())
    ax1.set_ylim(-25, 25)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_xlim(time.min(),time.max()+0.5)
    ax1.set_ylabel("Velocity (m/s)")
    ax2.plot(time, vcomp, 'b-')
    ax2.plot(time, ucomp, 'g-')
    ax2.set_xlim(time.min(),time.max()+0.5)
    ax2.set_ylim(-25, 25)
    ax2.set_xlabel("Date (UTC)")
    ax2.set_ylabel("Velocity (m/s)")
    ax2.xaxis.set_major_locator(MonthLocator())
    ax2.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax2.xaxis.set_minor_locator(DayLocator())
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('both')

    #fig1.autofmt_xdate()    
    # Set legend location - See: http://matplotlib.org/users/legend_guide.html#legend-location
    leg2 = plt.legend(['v','u'],loc='upper left')
    leg2._drawFrame=False
    DefaultSize = fig1.get_size_inches()
    fig1.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )

    fig1.suptitle("Ave Wind data for:   " + data_source, fontsize=12)
    # Save figure (without 'white' borders)
    plt.savefig('images/'+ station_name.lower() +'_' + data_source + '_timeseries.png', bbox_inches='tight', dpi = (100))
    plt.close(fig1) 
    
def dir_mag_hist(data1, data2, data3):

        #set mag bounds
    upper_lim = [2,8,100]
    lower_lim = [0,2,8]
    
    fig = plt.figure()
    ax = plt.subplot(311)
    
    plt.hist([data1[1][(data1[0] > lower_lim[0]) & (data1[0] <= upper_lim[0])], \
    data2[1][(data2[0] > lower_lim[0]) & (data2[0] <= upper_lim[0])], \
    data3[1][(data3[0] > lower_lim[0]) & (data3[0] <= upper_lim[0])]], np.arange(0,405,45), color=[(0,0,0),(.5,.5,.5),(1,1,1)])
    plt.show()

    ax = plt.subplot(312)
    
    plt.hist([data1[1][(data1[0] > lower_lim[1]) & (data1[0] <= upper_lim[1])], \
    data2[1][(data2[0] > lower_lim[1]) & (data2[0] <= upper_lim[1])], \
    data3[1][(data3[0] > lower_lim[1]) & (data3[0] <= upper_lim[1])]], np.arange(0,405,45), color=[(0,0,0),(.5,.5,.5),(1,1,1)])
    plt.show()
    
    ax = plt.subplot(313)
    
    plt.hist([data1[1][(data1[0] > lower_lim[2]) & (data1[0] <= upper_lim[2])], \
    data2[1][(data2[0] > lower_lim[2]) & (data2[0] <= upper_lim[2])], \
    data3[1][(data3[0] > lower_lim[2]) & (data3[0] <= upper_lim[2])]], np.arange(0,405,45), color=[(0,0,0),(.5,.5,.5),(1,1,1)])
    plt.show()    

    
"""---------------------------- Main Routine-------------------------------------------"""

"""------Ingest 1D Data--------"""

### NARR Data has the following boundary corners:
# Lambert Conformal
# 12.2N;133.5W, 54.5N; 152.9W, 57.3N; 49.4W ,14.3N;65.1W

year_long = '1997'
year_short = year_long[2:]
moor_sta_long = 'inner1'
moor_sta_short = 'ifm1'

NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/3hourly/'
NCEP = '/Users/bell/Data_Local/Reanalysis_Files/NCEPV2/6hourly/'
infile_narr = [NARR + 'uwnd.10m.'+year_long+'.nc', NARR + 'vwnd.10m.'+year_long+'.nc']
infile_ncep = [NCEP + 'uwnd.10m.gauss.'+year_long+'.nc', NCEP + 'vwnd.10m.gauss.'+year_long+'.nc']

### Grab grid points for future slicing - assume grid is same in all model output
narrlat_lon = latlon_grid(infile_narr[0])
nceplat_lon = latlon_grid(infile_ncep[0])

multifile=False
if multifile:
    try: #two possible file types amp and wpak
        MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/' + moor_sta_long.lower() + \
                        '/' + year_short + moor_sta_short +'*_amp.nc'
        MooringMetData, Mooring_params = from_netcdf_mf(MooringFile)
        MooringMetData = missing2nan(MooringMetData, Mooring_params)
    except RuntimeError:
        MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/' + moor_sta_long.lower() + \
                        '/' + year_short + moor_sta_short +'*_wpak.nc'
        MooringMetData, Mooring_params = from_netcdf_mf(MooringFile)
        MooringMetData = missing2nan(MooringMetData, Mooring_params)
else:
    try: #two possible file types amp and wpak
        MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/' + moor_sta_long.lower() + \
                    '/' + year_short + moor_sta_short +'a_amp.nc'
        MooringMetData, Mooring_params = from_netcdf(MooringFile)
        MooringMetData = missing2nan(MooringMetData, Mooring_params)
    except RuntimeError:
        MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/' + moor_sta_long.lower() + \
                    '/' + year_short + moor_sta_short +'a_wpak.nc'
        MooringMetData, Mooring_params = from_netcdf(MooringFile)
        MooringMetData = missing2nan(MooringMetData, Mooring_params)

MooringTime = date2pydate(MooringMetData['time'], MooringMetData['time2'], file_flag='EPIC')
sta_lat = MooringMetData['latitude'][0]
sta_long = MooringMetData['longitude'][0]

### force comparison location
#sta_lat = 56.
#sta_long = 165.

#Find NCEP and NARR nearest point to mooring
narrpt = sphered.nearest_point([sta_lat,-1 * sta_long],narrlat_lon['lat'],narrlat_lon['lon'], '2d')
nceppt = sphered.nearest_point([sta_lat,-1 * sta_long],nceplat_lon['lat'],nceplat_lon['lon']-360., '1d') #grid shift too

#Read in NARR and NCEP data for location chosen
NARR_uwind = from_netcdf_1dsplice(infile_narr[0], None, narrpt[3], narrpt[4])
NARR_vwind = from_netcdf_1dsplice(infile_narr[1], None, narrpt[3], narrpt[4])
NARRTime = date2pydate(NARR_uwind['time'], file_flag='NARR')

NCEP_uwind = from_netcdf_1dsplice(infile_ncep[0], 0, nceppt[3], nceppt[4])
NCEP_vwind = from_netcdf_1dsplice(infile_ncep[1], 0, nceppt[3], nceppt[4])
NCEPTime = date2pydate(NCEP_uwind['time'], file_flag='NCEP')

### calculate 6hr averages for all datasets using NARR time base
time_bin = 6.
time_str = str(time_bin) + 'hr'
NARRDaily_uwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), NARRTime, NARR_uwind['uwnd'], time_base=time_bin)
NARRDaily_vwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), NARRTime, NARR_vwind['vwnd'], time_base=time_bin)
NCEPDaily_uwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), NCEPTime, NCEP_uwind['uwnd'], time_base=time_bin)
NCEPDaily_vwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), NCEPTime, NCEP_vwind['vwnd'], time_base=time_bin)

MooringDaily_uwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), MooringTime, \
    wind_power_law(MooringMetData['WU_422'][:,0,0,0], correction=True), time_base=time_bin)
MooringDaily_vwnd = hourly_2_ave(NARRTime.min(),NARRTime.max(), MooringTime, \
    wind_power_law(MooringMetData['WV_423'][:,0,0,0], correction=True), time_base=time_bin)

"""---------------------------- Data Manipulation Routines-----------------------------"""

NARR_wind_mag = np.sqrt(NARRDaily_vwnd['mean']**2. + NARRDaily_uwnd['mean']**2.)
NARR_wind_dir_math = np.rad2deg(np.arctan2(NARRDaily_vwnd['mean'] , NARRDaily_uwnd['mean'])) 

NCEP_wind_mag = np.sqrt(NCEPDaily_vwnd['mean']**2. + NCEPDaily_uwnd['mean']**2.)
NCEP_wind_dir_math = np.rad2deg(np.arctan2(NCEPDaily_vwnd['mean'] , NCEPDaily_uwnd['mean'])) 

Mooring_wind_mag = np.sqrt(MooringDaily_uwnd['mean']**2. + MooringDaily_vwnd['mean']**2.)
Mooring_wind_dir_math = np.rad2deg(np.arctan2(MooringDaily_vwnd['mean'] , MooringDaily_uwnd['mean'])) 

# mask when mooring wasn't available
t_ind = ~np.isnan(Mooring_wind_mag) 

### Calculate +-flow and x-flow rotating along coast (~0 degrees bearing for no coordinate shift)
(NARRalong, NARRcross) = rotate_coord(0., NARR_wind_mag, NARR_wind_dir_math)
(NCEPalong, NCEPcross) = rotate_coord(0., NCEP_wind_mag, NCEP_wind_dir_math)
(MOORalong, MOORcross) = rotate_coord(0., Mooring_wind_mag, Mooring_wind_dir_math)

### normalize data
normalize = False
if normalize:
    NARRalong = (NARRalong ) / (NARR_wind_mag)
    NCEPalong = (NCEPalong ) / (NCEP_wind_mag)
    MOORalong = (MOORalong ) / (Mooring_wind_mag)
    NARRcross = (NARRcross ) / (NARR_wind_mag)
    NCEPcross = (NCEPcross ) / (NCEP_wind_mag)
    MOORcross = (MOORcross ) / (Mooring_wind_mag)
    
"""---------------------------- Plotting Routines--------------------------------------"""

### standard wind / time plots
#   NARR
quiver_timeseries(NARRDaily_uwnd['time'],NARRDaily_uwnd['mean'],NARRDaily_vwnd['mean'],NARR_wind_mag,'NARR', moor_sta_long.lower() )
quiver_timeseries(NCEPDaily_uwnd['time'],NCEPDaily_uwnd['mean'],NCEPDaily_vwnd['mean'],NCEP_wind_mag,'NCEP', moor_sta_long.lower() )
quiver_timeseries(MooringDaily_uwnd['time'],MooringDaily_uwnd['mean'],MooringDaily_vwnd['mean'],Mooring_wind_mag,'' + moor_sta_long.lower() + '', moor_sta_long.lower() )

###split timeseries
# 56 represents 2 weeks for 12 hour ave data
len2split = 56
len2split = int(len( MOORalong[t_ind] ))
split_int = np.ceil( len( MOORalong[t_ind] ) / len2split )
print split_int

for add_ind in range(0,int(split_int)):
    split_ind = range(add_ind*len2split,(add_ind*len2split)+len2split)
    if (split_ind[-1] > len( MOORalong[t_ind] )):
        split_ind = range(add_ind*len2split,len( MOORalong[t_ind] ))
        print "The last group has the following number of datapoints: " + np.str( len(split_ind) )
    print "Group " + np.str(add_ind)

    """ Most relevant plots below... along/U-comp coorelations"""
    ### Along/Cross Shore comparisons Mooring vs NARR/NCEP
    #  for entire year (mark mooring specific times)
    fig = plt.figure(6)
    #text locations
    right = 0.05
    top = .95

    (slope, intercept, r_value, p_value, std_err) = lin_fit(MOORalong[t_ind][split_ind], NARRalong[t_ind][split_ind])
    print "Regression stats for V-comp (along) Mooring v NARR are: %s %s " % (slope, r_value**2)
    (coor_mag, coor_angle)  = comp_corr((MOORcross[t_ind][split_ind],MOORalong[t_ind][split_ind]),(NARRcross[t_ind][split_ind],NARRalong[t_ind][split_ind]))
    print "NARR Complex correlation mag - %s and dir - %s" % (coor_mag, coor_angle)

    ax = plt.subplot(221)
    p1 = ax.plot(MOORalong[t_ind][split_ind], NARRalong[t_ind][split_ind])
    plt.setp(p1,color='r', marker='+', linestyle='')
    p2 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
    plt.setp(p2,'color','k','linestyle','--')
    p3 = ax.plot(np.arange(-15,16,5),(slope * np.arange(-15,16,5) + intercept) )
    plt.setp(p3,'color','k','linestyle','-.')
    ax.text(right, top, r"${r^2}$: %0.2f" % (r_value**2.),
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes, size=10)
    ax.set_xticks(np.arange(-15,16,5))
    ax.set_yticks(np.arange(-15,16,5))
    ax.set_xlim((-15,15))
    ax.set_ylim((-15,15))
    ax.set_xlabel(time_str + ' ' + moor_sta_long.lower() + ' V-comp Flow (m/s)')
    ax.set_ylabel(time_str + ' NARR V-comp Flow (m/s)')

    (slope, intercept, r_value, p_value, std_err) = lin_fit(MOORalong[t_ind][split_ind], NCEPalong[t_ind][split_ind])
    print "Regression stats for V-comp (along) Mooring v NCEP are: %s %s " % (slope, r_value**2)
    (coor_mag, coor_angle)  = comp_corr((MOORcross[t_ind][split_ind],MOORalong[t_ind][split_ind]),(NCEPcross[t_ind][split_ind],NCEPalong[t_ind][split_ind]))
    print "NCEP Complex correlation mag - %s and dir - %s" % (coor_mag, coor_angle)

    ax = plt.subplot(223)
    p1 = ax.plot(MOORalong[t_ind][split_ind], NCEPalong[t_ind][split_ind])
    plt.setp(p1,color='r', marker='+', linestyle='')
    p2 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
    plt.setp(p2,'color','k','linestyle','--')

    p3 = ax.plot(np.arange(-15,16,5),(slope * np.arange(-15,16,5) + intercept) )
    plt.setp(p3,'color','k','linestyle','-.')
    ax.text(right, top, r"${r^2}$: %0.2f" % (r_value**2.),
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes, size=10)

    ax.set_yticks(np.arange(-15,16,5))
    ax.set_xticks(np.arange(-15,16,5))
    ax.set_xlim((-15,15))
    ax.set_ylim((-15,15))
    ax.set_xlabel(time_str + ' ' + moor_sta_long.lower() + ' V-comp Flow (m/s)')
    ax.set_ylabel(time_str + ' NCEP V-comp Flow (m/s)')

    (slope, intercept, r_value, p_value, std_err) = lin_fit(MOORcross[t_ind][split_ind], NARRcross[t_ind][split_ind])
    print "Regression stats for U-comp Mooring (across) v NARR are: %s %s " % (slope, r_value**2)

    ax = plt.subplot(222)
    p1 = ax.plot(MOORcross[t_ind][split_ind], NARRcross[t_ind][split_ind])
    plt.setp(p1,color='r', marker='+', linestyle='')
    p2 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
    plt.setp(p2,'color','k','linestyle','--')

    p3 = ax.plot(np.arange(-15,16,5),(slope * np.arange(-15,16,5) + intercept) )
    plt.setp(p3,'color','k','linestyle','-.')
    ax.text(right, top, r"${r^2}$: %0.2f" % (r_value**2.),
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes, size=10)

    ax.set_xticks(np.arange(-15,16,5))
    ax.set_yticks(np.arange(-15,16,5))
    ax.set_xlim((-15,15))
    ax.set_ylim((-15,15))
    ax.set_xlabel(time_str + ' ' + moor_sta_long.lower() + ' U-comp Flow (m/s)')
    ax.set_ylabel(time_str + ' NARR U-comp Flow (m/s)')

    (slope, intercept, r_value, p_value, std_err) = lin_fit(MOORcross[t_ind][split_ind], NCEPcross[t_ind][split_ind])
    print "Regression stats for U-comp Mooring (across) v NCEP are: %s %s " % (slope, r_value**2)

    ax = plt.subplot(224)
    p1 = ax.plot(MOORcross[t_ind][split_ind], NCEPcross[t_ind][split_ind])
    plt.setp(p1,color='r', marker='+', linestyle='')
    p2 = ax.plot(np.arange(-15,16,5),np.arange(-15,16,5))
    plt.setp(p2,'color','k','linestyle','--')

    p3 = ax.plot(np.arange(-15,16,5),(slope * np.arange(-15,16,5) + intercept) )
    plt.setp(p3,'color','k','linestyle','-.')
    ax.text(right, top, r"${r^2}$: %0.2f" % (r_value**2.),
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes, size=10)

    ax.set_xticks(np.arange(-15,16,5))
    ax.set_yticks(np.arange(-15,16,5))
    ax.set_xlim((-15,15))
    ax.set_ylim((-15,15))
    ax.set_xlabel(time_str + ' ' + moor_sta_long.lower() + ' U-comp Flow (m/s)')
    ax.set_ylabel(time_str + ' NCEP U-comp Flow (m/s)')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )
    plt.savefig('images/' + moor_sta_long.lower() + '_alongacross_comp.png', bbox_inches='tight', dpi = (100))
    plt.close()


### Plot geolocations of datasets
plot_loc = True
if plot_loc:

    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=46, \
        urcrnrlat=62,llcrnrlon=-180,urcrnrlon=-140, lat_ts=45)
    
    lon_ncep, lat_ncep = np.meshgrid(-1. * nceplat_lon['lon'],nceplat_lon['lat'])
    x, y = m(lon_ncep, lat_ncep)

    #NARR - array given to define grid
    x_narr, y_narr = m(narrlat_lon['lon'],narrlat_lon['lat'])

    # Mooring Data
    x_moor, y_moor = m(-1. * sta_long,sta_lat)

    #ETOPO 5 contour data 
    ex, ey = m(elons, elats)
    CS = m.contourf(ex,ey,topoin, levels=range(250,5000,250), cmap='gray_r', alpha=.75) #colors='black'
    CS = m.contour(ex,ey,topoin, levels=range(250,5000,250), linewidths=0.2, colors='black', alpha=.75) #
    CS = m.contour(ex,ey,topoin, levels=[-1000, -200, -100], linestyle='--', linewidths=0.2, colors='black', alpha=.75) #
    plt.clabel(CS, inline=1, fontsize=8, fmt='%1.0f')
    
    #plot points
    m.scatter(x,y,20,marker='+',color='r', alpha=.75)
    m.scatter(x_narr,y_narr,20,marker='x',color='b', alpha=.75)
    m.scatter(x_moor,y_moor,20,marker='o',color='g', alpha=.75)

    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(46,66,4.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-180,-140,5.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )

    plt.savefig('images/' + moor_sta_long+ '_region.png', bbox_inches='tight', dpi = (100))
    plt.close()

