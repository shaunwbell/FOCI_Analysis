#!/usr/bin/env

"""
Program: NARR_anom.py

Purpose: create anomaly data from montly means and monthly longterm means 

DataSource: NARR data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/ 

Input/Arguments: Latitude and Longitude, raw data file
 
"""
#System Stack
import datetime
import sys

#Science Stack
from netCDF4 import Dataset
import numpy as np

# User Stack
import utilities.ncutilities as ncutil
import general_utilities.haversine as sphered

#Visual Packages
import matplotlib as mpl
#mpl.use('Agg') 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 03, 25)
__modified__ = datetime.datetime(2014, 03, 25)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'M2', 'SST','Bering'

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

def ncreadfile_dic_oneday(nchandle, params,doy,coords):
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
            if v in coords:
                data[v] = nchandle.variables[v][:]
            else:
                data[v] = nchandle.variables[v][doy]

        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)

def latlon_grid(infile):
    nchandle = ncutil.ncopen(infile)
    lat_lon = get_geocoords(nchandle)
    ncutil.ncclose(nchandle)

    return (lat_lon)

def get_geocoords(nchandle, lat='lat', lon='lon'):
    
    data = {}
    
    for j, v in enumerate([lat, lon]):
        data[v] = nchandle.variables[v][:]
        
    return (data)
    
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
        

    elif file_flag == 'NARR':
        """ Hours since 1800-1-1"""
        base_date=datetime.datetime.strptime('1800-01-01','%Y-%m-%d').toordinal()
        python_time = file_time / 24. + base_date
    elif file_flag == 'NCEP_days':
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

"""------------------------------------- MAPS -----------------------------------------"""

def etopo5_data():
    """ read in etopo5 topography/bathymetry. """
    file = 'data/etopo5.nc'
    etopodata = Dataset(file)
    
    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()
    
    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0
    
    #lons, lats = np.meshgrid(lons, lats)
    
    return(topoin, lats, lons)

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx

"""------------------------------- MAIN--------------------------------------------"""

anomperiod='daily'
ncvar = 'air'
ncfile_ltmean = '/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'+ ncvar +'.2m.day.ltm.nc'
ncfile_mean = '/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'+ ncvar +'.2m.2015.nc'

nchandle = Dataset(ncfile_ltmean,'r')
params = get_vars(nchandle)

data_lt = ncreadfile_dic(nchandle,params)
time_lt = date2pydate(data_lt['time'], file_flag='NARR')
nchandle.close()

nchandle = Dataset(ncfile_mean,'r')
params = get_vars(nchandle)

data = ncreadfile_dic(nchandle,params)
time = date2pydate(data['time'], file_flag='NARR')
time_label = [datetime.datetime.fromordinal(int(x)).strftime("%Y-%m-%d") for x in time]
nchandle.close()

### Grab grid points for future slicing - assume grid is same in all model output
lat_lon = latlon_grid(ncfile_mean)


station_name = ['Between M2-M4']
sta_lat = [57.31867]
sta_long = [166.3232]

'''
station_name = ['UnimakPass - Shumigan Downstream']
sta_lat = [54.5]
sta_long = [161.0]
'''
#Find NARR nearest point to moorings - haversine formula
# NARR data is -180->180 (positive east), Moorings are usually expressed +W for FOCI
station_pt = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '2d')
station_modelpt = [lat_lon['lat'][station_pt[3],station_pt[4]],lat_lon['lon'][station_pt[3],station_pt[4]]]

print "station nearest point to {0}, {1} which is lat:{2} , lon:{3}".format(sta_lat[0],sta_long[0], station_modelpt[0], station_modelpt[1])
    
"""------------------------------------------------------------------------------------"""
### Determine anomalies by finding the appropriate month to subtract off from the longterm data
#       anom_monthly = mon.mean - mon.ltm

anom = []
if anomperiod=='monthly':
    for index, each_date in enumerate(time_label):
        month_ind = int(each_date.split('-')[1]) - 1

        print "{0},{1},{2},{3}".format(each_date, data[ncvar][index,station_pt[3],station_pt[4]], \
            data_lt[ncvar][month_ind,station_pt[3],station_pt[4]], (data[ncvar][index] - data_lt[ncvar][month_ind])[station_pt[3],station_pt[4]])
    
        anom = np.hstack((anom, (data[ncvar][index] - data_lt[ncvar][month_ind])[station_pt[3],station_pt[4]]) )
elif anomperiod=='daily':
    for index, each_date in enumerate(time_label):

        print "{0},{1},{2},{3}".format(each_date, data[ncvar][index,station_pt[3],station_pt[4]], \
            data_lt[ncvar][index,station_pt[3],station_pt[4]], (data[ncvar][index] - data_lt[ncvar][index])[station_pt[3],station_pt[4]])
    
        anom = np.hstack((anom, (data[ncvar][index] - data_lt[ncvar][index])[station_pt[3],station_pt[4]]) )
else:
    print "No valid time base for data"

   
plot_geoloc = False
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
