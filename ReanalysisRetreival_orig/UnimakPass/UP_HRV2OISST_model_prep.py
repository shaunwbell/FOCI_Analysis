#!/usr/bin/env

"""
 UP_HRV2OISST_model_prep.py
 

 DataSource: ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/

 NOAA High Resolution SST data provided by the NOAA/OAR/ESRL PSD, 
  Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/ 
"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num

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
        ncinstance = ncutil.EPIC_NC_SST(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time[0]))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time1=time[0], time2=time[1], latitude=lat_lon[0], longitude=-1 * lat_lon[1], \
            depth_level=0. )
        ncinstance.add_data('T_25', data[0])
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
        """ days since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time + base_date
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
NARR = '/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/'
infile = [NARR + 'sst.day.anom.1981.v2.nc']


### Grab grid points for future slicing - assume grid is same in all model output
lat_lon = latlon_grid(infile[0])

#stn    ['1','2']
station_name = ['ShumaginDown']
sta_lat = [54.5]
sta_long = [161]

#Find NCEP nearest point to moorings - haversine formula
# NCEP data is 0->360 (positive east), Moorings are usually expressed +W for FOCI
stn1_pt = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '1d')
stn1_modelpt = [lat_lon['lat'][stn1_pt[3]],lat_lon['lon'][stn1_pt[4]]]

print "stn1 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[0], sta_long[0], stn1_modelpt[0], stn1_modelpt[1])

stn1_modelpt[1] = -1.*((180 - stn1_modelpt[1]) + 180)
print "thus converting lon to degrees W positive {0}".format(stn1_modelpt[1])

#loop over all requested data   
years = range(1981,2017)

for yy in years:
    # retrieve only these location's data
    # uwnd
    infile = NARR + 'sst.day.anom.'+ str(yy) + '.v2.nc'
    print "Working on file " + infile
    stn1_data = from_netcdf_1dsplice(infile, None, stn1_pt[3], stn1_pt[4])

    stn1_sst = stn1_data['anom'] 
    

    #convert to EPIC time
    pydate = date2pydate(stn1_data['time'], file_flag='NARR')
    epic_time, epic_time1 = pydate2EPIC(pydate)

    # output u,v wind components from model grid points
    save_to_nc = True
    if save_to_nc:
        # write to NetCDF
        outfile = 'data/NOAA_OI_SST_V2_anom_stn1_' + str(yy) + '.nc'
        print "Writing to Epic NetCDF " + outfile
        write2epic( outfile, station_name[0], [epic_time, epic_time1], stn1_modelpt, [stn1_sst,])
    
    
plot_geoloc = True
if plot_geoloc:
    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=52, \
        urcrnrlat=58,llcrnrlon=-165,urcrnrlon=-155, lat_ts=45)

    # Mooring Data
    x_moor, y_moor = m(-1. * sta_long[0],sta_lat[0])
    x_close, y_close = m(stn1_modelpt[1], stn1_modelpt[0])

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


    


