#!/usr/bin/env

"""
 UP_Winds_NARR_model_prep.py
 
 Retrieve NARR winds for two locations:
    Upstream Shumagin Islands - 54.5N 161W
    Downstream Shumagin Islands - 55.5N 158.5W

 Filter NARR winds with a triangular filter (1/4, 1/2, 1/4) and output every 3hrs
 Provide U, V (calculations for along,across are available but not archived in ncfiles)

 Save in EPIC NetCDF standard
 
 #Added subsample routine that uses mod function (%) to pull out times that is divisible by xxx
"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np
from netCDF4 import Dataset

# User Stack
import general_utilities.haversine as sphered
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR','Unimak', 'Shumagin','3hr filtered', 'U,V','Winds', 'Gulf of Alaska'

"""------------------------General   Modules-------------------------------------------"""

def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
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
    
def latlon_grid(infile):
    nchandle = ncutil.ncopen(infile)
    lat_lon = ncutil.get_geocoords(nchandle)
    ncutil.ncclose(nchandle)

    return (lat_lon)

def write2epic( file_name, stationid, time, lat_lon, data ):
        ncinstance = ncutil.EPIC_NC(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time[0]))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time1=time[0], time2=time[1], latitude=lat_lon[0], longitude=-1 * lat_lon[1], \
            depth_level=10. )
        ncinstance.add_data('WU_422', data[0])
        ncinstance.add_data('WV_423', data[1])
        ncinstance.close()

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
    else:
        print "time flag not recognized"
        sys.exit()
        
    return np.array(python_time)
    
def pydate2EPIC(file_time):
    ref_time_py = datetime.datetime.toordinal(datetime.datetime(1968, 5, 23))
    ref_time_epic = 2440000

    offset = ref_time_epic - ref_time_py

    time1 = np.floor(file_time) + offset #truncate to get day and add 2440000 for true julian day
    
    time2 = ( file_time - np.floor(file_time) ) * (1000. * 60. * 60.* 24.) #milliseconds since 0000GMT
    
    return(time1, time2)

def pythondate2str(pdate):
    (year,month,day) = datetime.datetime.fromordinal(int(pdate)).strftime('%Y-%b-%d').split('-')
    delta_t = pdate - int(pdate)
    dhour = str(int(np.floor(24 * (delta_t))))
    dmin = str(int(np.floor(60 * ((24 * (delta_t)) - np.floor(24 * (delta_t))))))
    dsec = str(int(np.floor(60 * ((60 * ((24 * (delta_t)) - np.floor(24 * (delta_t)))) - \
                    np.floor(60 * ((24 * (delta_t)) - np.floor(24 * (delta_t))))))))
                    
    #add zeros to time
    if len(dhour) == 1:
        dhour = '0' + dhour
    if len(dmin) == 1:
        dmin = '0' + dmin
    if len(dsec) == 1:
        dsec = '0' + dsec
                

    return year + '-' + month + '-' + day + ' ' + dhour+':'+dmin+':'+dsec 
    
"---"

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
    
def triangle_smoothing(data_in):
    weights=np.array([0.25,0.5,0.25])
    filtered_data = np.convolve(data_in,np.array(weights),'same') #edge effects
    
    return filtered_data

"""------------------------- Topo   Modules -------------------------------------------"""

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
        
"""------------------------- Main   Modules -------------------------------------------"""

### list of files
NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/3hourly/'
infile = [NARR + 'uwnd.10m.2003.nc']


### Grab grid points for future slicing - assume grid is same in all model output
lat_lon = latlon_grid(infile[0])

#stn    ['1','2']
station_name = ['ShumaginUp','ShumaginDown']
sta_lat = [54.75,55.5]
sta_long = [161,158.5]

#Find NARR nearest point to moorings - haversine formula
# NARR data is -180->180 (positive east), Moorings are usually expressed +W for FOCI
stn1_pt = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '2d')
stn2_pt = sphered.nearest_point([sta_lat[1],-1 * sta_long[1]],lat_lon['lat'],lat_lon['lon'], '2d')
stn1_modelpt = [lat_lon['lat'][stn1_pt[3],stn1_pt[4]],lat_lon['lon'][stn1_pt[3],stn1_pt[4]]]
stn2_modelpt = [lat_lon['lat'][stn2_pt[3],stn2_pt[4]],lat_lon['lon'][stn2_pt[3],stn2_pt[4]]]

print "stn1 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[0], sta_long[0], stn1_modelpt[0], stn1_modelpt[1])
    
print "stn2 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[1], sta_long[1], stn2_modelpt[0], stn2_modelpt[1])    

#loop over all requested data   
years = range(1980,2016)

for yy in years:
    # retrieve only these location's data
    # uwnd
    infile = NARR + 'uwnd.10m.'+ str(yy) + '.nc'
    print "Working on file " + infile
    stn1_data = from_netcdf_1dsplice(infile, None, stn1_pt[3], stn1_pt[4])
    stn2_data = from_netcdf_1dsplice(infile, None, stn2_pt[3], stn2_pt[4])

    #filter data
    stn1u_f = triangle_smoothing(stn1_data['uwnd'])
    stn2u_f = triangle_smoothing(stn2_data['uwnd'])

    stn1u = stn1_data['uwnd']
    stn2u = stn2_data['uwnd'] 
    
    # retrieve only these location's data
    # vwnd
    infile = NARR + 'vwnd.10m.'+ str(yy) + '.nc'
    print "Working on file " + infile
    stn1_data = from_netcdf_1dsplice(infile, None, stn1_pt[3], stn1_pt[4])
    stn2_data = from_netcdf_1dsplice(infile, None, stn2_pt[3], stn2_pt[4])

    #filter data
    stn1v_f = triangle_smoothing(stn1_data['vwnd'])
    stn2v_f = triangle_smoothing(stn2_data['vwnd'])

    stn1v = stn1_data['vwnd']
    stn2v = stn2_data['vwnd']

    #rotate to shore (Along/Across)
    #stn1
    '''
    NARR_wind_mag = np.sqrt(stn1u**2. + stn1v**2.)
    NARR_wind_dir_math = np.rad2deg(np.arctan2(stn1v, stn1u)) 
    (NARRalong, NARRcross) = rotate_coord(135., NARR_wind_mag, NARR_wind_dir_math)
    '''
    
    #convert to EPIC time
    pydate = date2pydate(stn1_data['time'], file_flag='NARR')
    epic_time, epic_time1 = pydate2EPIC(pydate)
    
    ###
    #output 0,6,12,18 UTC
    #subsample data
    time_ind = np.where(pydate%0.25 == 0)[0]
    
    # output u,v wind components from model grid points
    save_to_nc = True
    if save_to_nc:
        # write to NetCDF
        outfile = 'data/NARR_stn1_' + str(yy) + '.nc'
        print "Writing to Epic NetCDF " + outfile
        write2epic( outfile, station_name[1], [epic_time[time_ind], epic_time1[time_ind]], stn1_modelpt, [stn1u_f[time_ind], stn1v_f[time_ind]])
        outfile = 'data/NARR_stn2_' + str(yy) + '.nc'
        print "Writing to Epic NetCDF " + outfile
        write2epic( outfile, station_name[0], [epic_time[time_ind], epic_time1[time_ind]], stn2_modelpt, [stn2u_f[time_ind], stn2v_f[time_ind]])
    
    output2screen = False
    if output2screen:
        print"Date/Time, Across (m/s), Along(m/s)\n"
        for i,v in enumerate(pydate):
            print "{0}, {1}, {2}".format(pythondate2str(v), NARRcross[i],NARRalong[i])
    
plot_geoloc = True
if plot_geoloc:
    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=52, \
        urcrnrlat=58,llcrnrlon=-165,urcrnrlon=-155, lat_ts=45)

    # Mooring Data
    x_moor, y_moor = m([-1. * sta_long[0], -1. * sta_long[1]],sta_lat)
    x_close, y_close = m([stn1_modelpt[1],stn2_modelpt[1]], [stn1_modelpt[0],stn2_modelpt[0]])

    #ETOPO 5 contour data 
    ex, ey = m(elons, elats)
    CS = m.contourf(ex,ey,topoin, levels=range(250,5000,250), cmap='gray_r', alpha=.75) #colors='black'
    CS = m.contour(ex,ey,topoin, levels=range(250,5000,250), linewidths=0.2, colors='black', alpha=.75) #
    CS = m.contour(ex,ey,topoin, levels=[-1000, -200, -100], linestyle='--', linewidths=0.2, colors='black', alpha=.75) #
    
    #plot points
    m.scatter(x_close,y_close,20,marker='+',color='b')
    m.scatter(x_moor,y_moor,20,marker='o',color='g')

    m.drawcountries(linewidth=0.5)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(50,62,2.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-165,-145,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]) )

    plt.savefig('images/shumigans_region.png', bbox_inches='tight', dpi = (100))
    plt.close()


    


