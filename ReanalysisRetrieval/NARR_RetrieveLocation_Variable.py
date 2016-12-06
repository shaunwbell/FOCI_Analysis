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

#User Stack
from io_utils.EcoFOCI_netCDF_read import EcoFOCI_netCDF
from calc.EPIC2Datetime import EPIC2Datetime, get_UDUNITS
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

"""--------------------------------main Routines---------------------------------------"""


""" currently hard coded - variables and ranges """


### Grab grid points for future slicing - assume grid is same in all model output
NARR = '/Users/bell/Data_Local/Reanalysis_Files/NARR/3hourly/'
infile = [NARR + 'uwnd.10m.2010.nc']

lat_lon = latlon_grid(infile[0])

#stn    ['1','2']
station_name = ['C2']
sta_lat = [71.23]
sta_long = [164.221]

#Find NARR nearest point to moorings - haversine formula
# NARR data is -180->180 (positive east), Moorings are usually expressed +W for FOCI
station_1 = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '2d')
stn1_modelpt = [lat_lon['lat'][station_1[3],station_1[4]],lat_lon['lon'][station_1[3],station_1[4]]]

print "stn1 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[0], sta_long[0], stn1_modelpt[0], stn1_modelpt[1])
    

#loop over all requested data   
years = range(2010,2017)
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

"""-----------------------------------"""

#read in NARR datafile
ncfile = args.DataPath
df = EcoFOCI_netCDF(file)
global_atts = df.get_global_atts()
vars_dic = df.get_vars()
ncdata = df.ncreadfile_dic()
df.close()