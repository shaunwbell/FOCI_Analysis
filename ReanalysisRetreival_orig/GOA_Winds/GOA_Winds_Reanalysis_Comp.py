#!/usr/bin/env

"""
 GOA_Winds_Reanalysis_Comp.py
 
 Compare NARR Winds with NCEP V2 (with Mooring Winds)

 Using Anaconda packaged Python 
"""
#System Stack
import datetime

#Science Stack
import numpy as np
from netCDF4 import Dataset
from scipy import stats

# User Stack
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
__keywords__ = 'NCEP Reanalysis V2','NARR','GLOBEC3','Daily comparison', 'Winds', 'Gulf of Alaska'

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
        sys.exit()
        
    return np.array(python_time)
        
"""------------------------------- Stats/Math Modules --------------------------------------"""

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
    
def wind_power_law(comp_orig, height_obs=3., height_interp=10., correction=False):
    """simple power law wind adjustment
    default - 3m observations, 10m interpolated height"""
    
    if correction:
        wind_cor = comp_orig * (height_interp / height_obs)**(0.143)
    else:
        wind_cor = comp_orig
        
    return wind_cor
    

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
    ax2.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax2.xaxis.set_minor_locator(DayLocator())
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('both')

    fig1.autofmt_xdate()    
    # Set legend location - See: http://matplotlib.org/users/legend_guide.html#legend-location
    leg2 = plt.legend(['v','u'],loc='upper left')
    leg2._drawFrame=False
    DefaultSize = fig1.get_size_inches()
    fig1.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )

    fig1.suptitle("24hr ave Wind data for:   " + data_source, fontsize=12)
    # Save figure (without 'white' borders)
    plt.savefig('images/Globec_' + data_source + '_timeseries.png', bbox_inches='tight', dpi = (100))
    plt.close(fig1) 

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
        
"""---------------------------- Main Routine-------------------------------------------"""

"""------Ingest Data--------"""

year_long = '2001'
year_short = '01'

NCEPV2 = '/Users/bell/Data_Local/Reanalysis_Files/NCEPV2/daily_mean/'

NCEPV2_uwind, NCEPV2_uparams = from_netcdf(NCEPV2 + 'uwnd.10m.gauss.' + year_long + '.nc')
NCEPV2_vwind, NCEPV2_vparams = from_netcdf(NCEPV2 + 'vwnd.10m.gauss.' + year_long + '.nc')
NCEPTime = date2pydate(NCEPV2_uwind['time'], file_flag='NCEP')

NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'

NARR_uwind, NARR_uparams = from_netcdf(NARR + 'uwnd.10m.' + year_long + '.nc')
NARR_vwind, NARR_vparams = from_netcdf(NARR + 'vwnd.10m.' + year_long + '.nc')
NARRTime = date2pydate(NARR_uwind['time'], file_flag='NCEP')

### NARR Data has the following boundary corners:
# 12.2N;133.5W, 54.5N; 152.9W, 57.3N; 49.4W ,14.3N;65.1W
# Lambert Conformal

multifile=True
if multifile:
    MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/globec3/' + year_short + 'gbm3*_wpak.nc'
    MooringMetData, Mooring_params = from_netcdf_mf(MooringFile)
else:
    MooringFile = '/Users/bell/Data_Local/FOCI/Mooring/' + year_long + '/globec3/' + year_short + 'gbm3a_wpak.nc'
    MooringMetData, Mooring_params = from_netcdf(MooringFile)




MooringTime = date2pydate(MooringMetData['time'], MooringMetData['time2'], file_flag='EPIC')
# Use vectors, not deg/speed values for averages
MooringDaily_uwnd = hourly2daily(NARRTime.min(),NARRTime.max(), MooringTime, \
    wind_power_law(MooringMetData['WU_422'], correction=True))
MooringDaily_vwnd = hourly2daily(NARRTime.min(),NARRTime.max(), MooringTime, \
    wind_power_law(MooringMetData['WV_423'], correction=True))
sta_lat = MooringMetData['latitude'][0]
sta_long = MooringMetData['longitude'][0]

#-----> user set to force mooring location instead of using built in location (useful if you want
# to specify lat/lon for model comparison purposes
#sta_lat = 56.
#sta_long = 146.
"""---------------------------- Data Manipulation Routines-----------------------------"""

#NCEP V2 - vectors given to define grid
#shift grid from 0->360 to -360->0
#only the longitude is relevant, the index in the array remains the same
NCEP_uwnds,lons_ncep = shiftgrid(0.,NCEPV2_uwind['uwnd'],NCEPV2_uwind['lon'],start=False) 

#Find NCEP and NARR nearest point to mooring
nceppt = sphered.nearest_point([sta_lat,-1 * sta_long],NCEPV2_uwind['lat'],lons_ncep, '1d')
narrpt = sphered.nearest_point([sta_lat,-1 * sta_long],NARR_uwind['lat'],NARR_uwind['lon'], '2d')

### grid specific
# get met mag and dir (N -> 90deg and mag is with direction of flow ... mass/oceanographic standard

#-----> user set to turn on/off two point interpolation
twopt_ave = False
if twopt_ave: #average two points in latitude - closest to station and one to the immediate south
    NCEP_uwind = NCEPV2_uwind['uwnd'][:,:,nceppt[3]:nceppt[3]+2,nceppt[4]].mean(axis=2)
    NCEP_vwind = NCEPV2_vwind['vwnd'][:,:,nceppt[3]:nceppt[3]+2,nceppt[4]].mean(axis=2)
else:
    NCEP_uwind = NCEPV2_uwind['uwnd'][:,:,nceppt[3],nceppt[4]]
    NCEP_vwind = NCEPV2_vwind['vwnd'][:,:,nceppt[3],nceppt[4]]

NCEP_wind_mag = np.sqrt(NCEP_uwind**2. + NCEP_vwind**2.)
NCEP_wind_dir_math = np.rad2deg(np.arctan2(NCEP_vwind , NCEP_uwind)) 
NCEP_wind_dir = cart2wind(NCEP_wind_dir_math)

try: #some data gives extra paramter and puts data in a dict structure... others does not
    NARR_u = NARR_uwind['uwnd'].data[:,narrpt[3],narrpt[4]]
    NARR_v = NARR_vwind['vwnd'].data[:,narrpt[3],narrpt[4]]
    
except TypeError: #no .data parameter
    NARR_u = NARR_uwind['uwnd'][:,narrpt[3],narrpt[4]]
    NARR_v = NARR_vwind['vwnd'][:,narrpt[3],narrpt[4]]
    
NARR_wind_mag = np.sqrt(NARR_v**2. + NARR_u**2.) # corrected NARR_u -> NARR_v 9/11
NARR_wind_dir_math = np.rad2deg(np.arctan2(NARR_v , NARR_u)) 
NARR_wind_dir = cart2wind(NARR_wind_dir_math)

Mooring_wind_mag = np.sqrt(MooringDaily_uwnd['daily_mean']**2. + MooringDaily_vwnd['daily_mean']**2.)
Mooring_wind_dir_math = np.rad2deg(np.arctan2(MooringDaily_vwnd['daily_mean'] , MooringDaily_uwnd['daily_mean'])) 
Mooring_wind_dir = cart2wind(Mooring_wind_dir_math)

# mask when mooring wasn't available
t_ind = ~np.isnan(Mooring_wind_mag) & (Mooring_wind_mag < 100)

### Calculate +-flow and x-flow rotating along coast (~43 degrees bearing near Globec3 )
(NCEPalong, NCEPcross) = rotate_coord(137., NCEP_wind_mag[:,0], NCEP_wind_dir_math[:,0])
(NARRalong, NARRcross) = rotate_coord(137., NARR_wind_mag, NARR_wind_dir_math)
(MOORalong, MOORcross) = rotate_coord(137., Mooring_wind_mag, Mooring_wind_dir_math)

"""---------------------------- Plotting Routines--------------------------------------"""

### standard wind / time plots
# NARR
quiver_timeseries(NARRTime,NARR_u,NARR_v,NARR_wind_mag,'NARR')
quiver_timeseries(NCEPTime,NCEP_uwind[:,0],NCEP_vwind[:,0],NCEP_wind_mag[:,0],'NCEP_V2')
quiver_timeseries(MooringDaily_uwnd['daily_time'],MooringDaily_uwnd['daily_mean'],MooringDaily_vwnd['daily_mean'],Mooring_wind_mag,'GLOBEC3')


"""
### one-to-one NARR/NCEP vs Mooring magnitude in met coords for mooring times only
fig = plt.figure(2)
ax = plt.subplot(121)
p1 = ax.plot(Mooring_wind_mag[t_ind], NCEP_wind_mag[t_ind,0], label="vs NCEP V2")
plt.setp(p1,color='r', marker='.', linestyle='')
p2 = ax.plot(Mooring_wind_mag[t_ind], NARR_wind_mag[t_ind], label="vs NARR")
plt.setp(p2,color='b', marker='+', linestyle='')
p3 = ax.plot(np.arange(0,16,1),np.arange(0,16,1))
plt.setp(p3,'color','k','linestyle','--')
ax.set_ylabel('Daily Wind Speed m/s Reanalysis')
ax.set_xlabel('Daily Wind Speed m/s Globec3 Mooring')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='lower right')
ax = plt.subplot(122)
plt.hist([NCEP_wind_mag[t_ind,0], NARR_wind_mag[t_ind], Mooring_wind_mag[t_ind]],\
    np.arange(0,16,1), stacked=False, orientation='horizontal', color=['r','b','k'], rwidth=0.9)
#plt.tight_layout()
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_NARRNCEP_wscomp.png', bbox_inches='tight', dpi = (100))
plt.close()

### one-to-one NARR/NCEP vs Mooring direction in met coords for mooring times only
fig = plt.figure(3)
ax3 = plt.subplot(121)
p1 = ax3.plot(Mooring_wind_dir[t_ind], NCEP_wind_dir[t_ind,0], label="vs NCEP V2")
plt.setp(p1,color='r', marker='.', linestyle='')
p2 = ax3.plot(Mooring_wind_dir[t_ind], NARR_wind_dir[t_ind], label="vs NARR")
plt.setp(p2,color='b', marker='+', linestyle='')
p3 = ax3.plot(np.arange(0,361,45),np.arange(0,361,45))
plt.setp(p3,'color','k','linestyle','--')
ax3.set_xticks(np.arange(0,361,45))
ax3.set_yticks(np.arange(0,361,45))
ax3.set_ylabel('Daily Wind Direction Reanalysis')
ax3.set_xlabel('Daily Wind Direction Globec3')
handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc='upper center')
ax = plt.subplot(122)
plt.hist([NCEP_wind_dir[t_ind,0], NARR_wind_dir[t_ind], Mooring_wind_dir[t_ind]], \
    np.arange(0,361,30), stacked=False, orientation='horizontal', color=['r','b','k'], rwidth=0.9)
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_NARRNCEP_degcomp.png', bbox_inches='tight', dpi = (100))
plt.close()

### one-to-one NARRvsNCEP for entire year (mark mooring specific times)
fig = plt.figure(4)
ax = plt.subplot(121)
p1 = ax.plot(NARR_wind_mag, NCEP_wind_mag[:,0])
plt.setp(p1,color='r', marker='+', linestyle='')
p2 = ax.plot(NARR_wind_mag[t_ind], NCEP_wind_mag[t_ind,0])
plt.setp(p2,color='k', marker='.', linestyle='')
p3 = ax.plot(np.arange(0,26,5),np.arange(0,26,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(0,26,5))
ax.set_yticks(np.arange(0,26,5))
ax.set_xlabel('Daily Wind Speed m/s NARR')
ax.set_ylabel('Daily Wind Speed m/s NCEP V2')
ax = plt.subplot(122)
p1 = ax.plot(NARR_wind_dir, NCEP_wind_dir[:,0])
plt.setp(p1,color='r', marker='.', linestyle='')
plt.setp(p1,color='r', marker='+', linestyle='')
p2 = ax.plot(NARR_wind_dir[t_ind], NCEP_wind_dir[t_ind,0])
plt.setp(p2,color='k', marker='.', linestyle='')
p3 = ax.plot(np.arange(0,361,45),np.arange(0,361,45))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(0,361,45))
ax.set_yticks(np.arange(0,361,45))
ax.set_xlabel('Daily Wind Dir NARR')
ax.set_ylabel('Daily Wind Dir NCEP V2')
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_NARRNCEP_magdeg.png', bbox_inches='tight', dpi = (100))
plt.close()

### one-to-one U,V wind NARRvsNCEP for entire year (mark mooring specific times)
fig = plt.figure(41)
ax = plt.subplot(111)
p1 = ax.plot(NARR_u, NCEP_uwind[:,0])
plt.setp(p1,color='r', marker='+', linestyle='')
p2 = ax.plot(NARR_v, NCEP_vwind[:,0])
plt.setp(p2,color='k', marker='.', linestyle='')
p3 = ax.plot(np.arange(-25,26,5),np.arange(-25,26,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(-25,26,5))
ax.set_yticks(np.arange(-25,26,5))
ax.set_xlabel('Daily Wind Components (m/s) NARR')
ax.set_ylabel('Daily Wind Components (m/s) NCEP V2')
ax.legend(['U comp','V comp'], loc='lower right')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_NARRNCEP_uv.png', bbox_inches='tight', dpi = (100))
plt.close()

### one-to-one U,V wind NARRvsNCEP for entire year (mark mooring specific times)
fig = plt.figure(42)
ax = plt.subplot(121)
p1 = ax.plot(MooringDaily_uwnd['daily_mean'][t_ind], NCEP_uwind[t_ind,0])
plt.setp(p1,color='r', marker='+', linestyle='')
p2 = ax.plot(MooringDaily_vwnd['daily_mean'][t_ind], NCEP_vwind[t_ind,0])
plt.setp(p2,color='k', marker='.', linestyle='')
p3 = ax.plot(np.arange(-25,26,5),np.arange(-25,26,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(-25,26,5))
ax.set_yticks(np.arange(-25,26,5))
ax.set_xlabel('Daily Wind Components (m/s) Globec3')
ax.set_ylabel('Daily Wind Components (m/s) NCEP V2')
ax.legend(['U comp','V comp'], loc='lower right')

ax = plt.subplot(122)
p1 = ax.plot(MooringDaily_uwnd['daily_mean'][t_ind], NARR_u[t_ind])
plt.setp(p1,color='r', marker='+', linestyle='')
p2 = ax.plot(MooringDaily_vwnd['daily_mean'][t_ind], NARR_v[t_ind])
plt.setp(p2,color='k', marker='.', linestyle='')
p3 = ax.plot(np.arange(-25,26,5),np.arange(-25,26,5))
plt.setp(p3,'color','k','linestyle','--')
ax.set_xticks(np.arange(-25,26,5))
ax.set_yticks(np.arange(-25,26,5))
ax.set_xlabel('Daily Wind Components (m/s) Globec3')
ax.set_ylabel('Daily Wind Components (m/s) NARR')
ax.legend(['U comp','V comp'], loc='lower right')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_GlobecNARRNCEP_uv.png', bbox_inches='tight', dpi = (100))
plt.close()

### Mooring vs NARR/NCEP wind offset met coords for mooring times only
fig = plt.figure(5)
ax = plt.subplot(111)
plt.plot(Mooring_wind_dir[t_ind], Mooring_wind_dir[t_ind] - NARR_wind_dir[t_ind],'r.')
plt.plot(Mooring_wind_dir[t_ind], Mooring_wind_dir[t_ind] - NCEP_wind_dir[t_ind,0],'b.')
ax.set_xticks(np.arange(0,361,45))
ax.set_yticks(np.arange(-360,361,45))
ax.set_xlabel('Mooring Wind Direction')
ax.set_ylabel('Reanalysis vs Mooring Diff (Direction)')
ax.grid(True)
ax.legend(['NARR','NCEP V2'], loc='lower right')
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig('images/Globec_winddif.png', bbox_inches='tight', dpi = (100))
plt.close()
"""

""" Most relevant plots below... along/across shore coorelations"""
### Along/Cross Shore comparisons Mooring vs NARR/NCEP
#  for entire year (mark mooring specific times)
fig = plt.figure(6)
#text locations
right = 0.05
top = .95

(slope, intercept, r_value, p_value, std_err) = lin_fit(MOORalong[t_ind], NARRalong[t_ind])
print "Regression stats for Along Shore Mooring v NARR are: %s %s %s %s %s" % (slope, intercept, r_value, p_value, std_err)
(coor_mag, coor_angle)  = comp_corr((MOORcross[t_ind],MOORalong[t_ind]),(NARRcross[t_ind],NARRalong[t_ind]))
print "Complex correlation mag - %s and dir - %s" % (coor_mag, coor_angle)

ax = plt.subplot(221)
p1 = ax.plot(MOORalong[t_ind], NARRalong[t_ind])
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
ax.set_xlabel('Daily Globec3 Along-shore Flow (m/s)')
ax.set_ylabel('Daily NARR Along-shore Flow (m/s)')

(slope, intercept, r_value, p_value, std_err) = lin_fit(MOORalong[t_ind], NCEPalong[t_ind])
print "Regression stats for Along Shore Mooring v NCEP are: %s %s %s %s %s" % (slope, intercept, r_value, p_value, std_err)
(coor_mag, coor_angle)  = comp_corr((MOORcross[t_ind],MOORalong[t_ind]),(NCEPcross[t_ind],NCEPalong[t_ind]))
print "Complex correlation mag - %s and dir - %s" % (coor_mag, coor_angle)

ax = plt.subplot(223)
p1 = ax.plot(MOORalong[t_ind], NCEPalong[t_ind])
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
ax.set_xlabel('Daily Globec3 Along-shore Flow (m/s)')
ax.set_ylabel('Daily NCEP Along-shore Flow (m/s)')

(slope, intercept, r_value, p_value, std_err) = lin_fit(MOORcross[t_ind], NARRcross[t_ind])
print "Regression stats for Across Shore Mooring v NARR are: %s %s %s %s %s" % (slope, intercept, r_value, p_value, std_err)

ax = plt.subplot(222)
p1 = ax.plot(MOORcross[t_ind], NARRcross[t_ind])
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
ax.set_xlabel('Daily Globec3 Across-shore Flow (m/s)')
ax.set_ylabel('Daily NARR Across-shore Flow (m/s)')

(slope, intercept, r_value, p_value, std_err) = lin_fit(MOORcross[t_ind], NCEPcross[t_ind])
print "Regression stats for Across Shore Mooring v NCEP are: %s %s %s %s %s" % (slope, intercept, r_value, p_value, std_err)

ax = plt.subplot(224)
p1 = ax.plot(MOORcross[t_ind], NCEPcross[t_ind])
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
ax.set_xlabel('Daily Globec3 Across-shore Flow (m/s)')
ax.set_ylabel('Daily NCEP Across-shore Flow (m/s)')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )
plt.savefig('images/Globec_alongacross_comp.svg', bbox_inches='tight', dpi = (100))
plt.close()


### Plot geolocations of datasets
plot_loc = True
if plot_loc:

    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=55, \
        urcrnrlat=62,llcrnrlon=-155,urcrnrlon=-145, lat_ts=45)
    
    lon_ncep, lat_ncep = np.meshgrid(lons_ncep, NCEPV2_uwind['lat'])
    x, y = m(lon_ncep, lat_ncep)

    #NARR - array given to define grid
    x_narr, y_narr = m(NARR_uwind['lon'], NARR_uwind['lat'])

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
    m.drawparallels(np.arange(55,62,2.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-155,-145,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*1.5) )

    plt.savefig('images/Globec_region.png', bbox_inches='tight', dpi = (100))
    plt.close()
