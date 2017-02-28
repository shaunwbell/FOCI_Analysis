#!/usr/bin/env python

"""
 Background:
 --------
 NARR_RetrieveLocation_Variable.py
 
 
 Purpose:
 --------
 Routines to retrieve, output NARR data from a single point over time to combine for analysis
 
 History:
 --------

 2016-09-20 : Bell - simplify existing multiple routines for various locations into one package

"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np
from netCDF4 import num2date

#User Stack
from io_utils.EcoFOCI_netCDF_read import EcoFOCI_netCDF
from calc.EPIC2Datetime import EPIC2Datetime, get_UDUNITS, Datetime2EPIC
import calc.haversine as sphered

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 9, 20)
__modified__ = datetime.datetime(2016, 9, 20)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR'

    
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

    
def from_netcdf_1dsplice(infile, height_ind, lat_ind, lon_ind):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    df = EcoFOCI_netCDF(infile)
    nchandle = df._getnchandle_()

    params = df.get_vars() #gets all of them

    print "Parameters available: " 
    #print params
    
    ncdata = ncreadfile_dic_slice(nchandle, params, height_ind=height_ind, lat_ind=lat_ind, lon_ind=lon_ind)
    df.close()
    
    return ncdata
    
def get_geocoords(infile, lat='lat', lon='lon'):

    df = EcoFOCI_netCDF(infile)
    nchandle = df._getnchandle_()

    data = {}
    
    for j, v in enumerate([lat, lon]):
        data[v] = nchandle.variables[v][:]

    df.close()
        
    return (data)

def ncreadfile_dic_slice(nchandle, params, height_ind=None, lat_ind=None, lon_ind=None):
    """returns slice of data for all times but for specified height/lat/lon indicies"""
    data = {}
    if height_ind == None:
        for j, v in enumerate(params): 
            try: #check for nc variable
                    data[v] = nchandle.variables[v][:,lat_ind,lon_ind]

            except ValueError: #if parameter is not of expected dimensions
                data[v] = nchandle.variables[v][:]
    else:
        for j, v in enumerate(params): 
            try: #check for nc variable
                    data[v] = nchandle.variables[v][:,:,lat_ind,lon_ind]

            except ValueError: #if parameter is not of expected dimensions
                data[v] = nchandle.variables[v][:]

    return data    
    
"""--------------------------------main Routines---------------------------------------"""


""" currently hard coded - variables and ranges """


### Grab grid points for future slicing - assume grid is same in all model output
NARR = '/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'
infile = [NARR + 'uwnd.10m.2016.nc']

lat_lon = get_geocoords(infile[0])

#stn    ['1','2']
station_name = ['UP stn_1']
sta_lat = [54.5]
sta_long = [161.0]

#Find NARR nearest point to moorings - haversine formula
# NARR data is -180->180 (positive east), Moorings are usually expressed +W for FOCI
station_1 = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '2d')
stn1_modelpt = [lat_lon['lat'][station_1[3],station_1[4]],lat_lon['lon'][station_1[3],station_1[4]]]

print "stn1 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[0], sta_long[0], stn1_modelpt[0], stn1_modelpt[1])
    
"""
#loop over all requested data   
years = range(2010,2017)
years = ['mon.mean']
for yy in years:
    # retrieve only these location's data
    # uwnd
    infile = NARR + 'uwnd.10m.'+ str(yy) + '.nc'
    print "Working on file " + infile
    stn1_data = from_netcdf_1dsplice(infile, None, station_1[3], station_1[4])

    #filter data
    stn1u_f = triangle_smoothing(stn1_data['uwnd'])

    stn1u = stn1_data['uwnd']
    
    # retrieve only these location's data
    # vwnd
    infile = NARR + 'vwnd.10m.'+ str(yy) + '.nc'
    print "Working on file " + infile
    stn1_data = from_netcdf_1dsplice(infile, None, station_1[3], station_1[4])

    #filter data
    stn1v_f = triangle_smoothing(stn1_data['vwnd'])

    stn1v = stn1_data['vwnd']

    
    #convert to EPIC time
    #epic_time, epic_time1 = Datetime2EPIC(num2date(stn1_data['time'], "hours since 1800-1-1 00:00:0.0"))
    Datetime2EPIC(num2date(x, "hours since 1800-1-1 00:00:0.0")) for x in stn1_data['time']
    ###
    #output 0,6,12,18 UTC
    #subsample data
#    time_ind = np.where(pydate%0.25 == 0)[0]
    
    # output u,v wind components from model grid points
    save_to_nc = False
    if save_to_nc:
        # write to NetCDF
        outfile = 'data/NARR_stn1_' + str(yy) + '.nc'
        print "Writing to Epic NetCDF " + outfile
#        write2epic( outfile, station_name[1], [epic_time[time_ind], epic_time1[time_ind]], stn1_modelpt, [stn1u_f[time_ind], stn1v_f[time_ind]])
        write2epic( outfile, station_name[1], [epic_time, epic_time1], stn1_modelpt, [stn1u, stn1v])
"""
"""-----------using xarray---------"""
import pandas as pd
import xarray as xa

#index = [station_1[3],station_1[4]]
index=[195,76]
ufilein='/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/uwnd.10m.2016.nc'
udata = xa.open_dataset(ufilein, decode_cf=False)
udata = xa.decode_cf(udata,mask_and_scale=False)
dum = udata.uwnd[:443,195,76].resample('D', udata.time, how='mean')
print dum.to_pandas().to_csv()

vfilein='/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/vwnd.10m.2016.nc'
vdata = xa.open_dataset(vfilein, decode_cf=False)
vdata = xa.decode_cf(vdata,mask_and_scale=False)
dvm = vdata.vwnd[:443,195,76].resample('D', vdata.time, how='mean')
print dvm.to_pandas().to_csv()
