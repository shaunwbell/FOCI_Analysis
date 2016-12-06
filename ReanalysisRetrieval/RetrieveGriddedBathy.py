#!/usr/bin/env

"""
 RetrieveGriddedBathy.py
 
 Purpose:
 --------
 Choose Dataset (Etopo5, ArdemV2, Etopo1) and output bathymetry for selected lat/lon pair

 DataOptions - Dataset Source
 LocationOptions - Nearest Point / Interpolation


 Caveats:
 --------

 History:
 --------

"""

#System Stack
import datetime
import sys
import argparse

#Science Stack
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid

#User Stack
import calc.haversine as sphered

"""------------------------- Topo   Modules -------------------------------------------"""

def etopo5_data(region=None):
    """ read in etopo5 topography/bathymetry. """

    file = 'data/etopo5.nc'
    etopodata = Dataset(file)

    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()


    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0

    lons, lats = np.meshgrid(lons, lats)

    if region == 'Chukchi':
        latbounds = [60,90]
        lonbounds = [-180, -120]

        latind = np.ma.masked_outside(lats,latbounds[0],latbounds[1])
        lonind = np.ma.masked_outside(lons,lonbounds[0],lonbounds[1])
        ind = (latind.mask==True) &(latind.mask == lonind.mask)
        mtopoin = np.ma.masked_array(topoin,ind)
        lons = lons[ind]
        lats = lats[ind]

        
        return(mtopoin, latind, lonind)

    else:
        return(topoin, lats, lons)

def ardemV2_data(region=None):
    """ read in ardemV2 topography/bathymetry. """
    file = 'data/ARDEMV2.0.nc'
    ardemv2data = Dataset(file)
    
    topoin = ardemv2data.variables['z'][:]
    lons = ardemv2data.variables['lon'][:]
    lats = ardemv2data.variables['lat'][:]
    ardemv2data.close()
    
    if region == 'Chukchi':
        latbounds = [60,90]
        lonbounds = [-180, -120]

        latind = np.ma.masked_outside(lats,latbounds[0],latbounds[1])
        lonind = np.ma.masked_outside(lons,lonbounds[0],lonbounds[1])
        ind = (latind.mask==True) &(latind.mask == lonind.mask)
        mtopoin = np.ma.masked_array(topoin,ind)
        lons = lons[ind]
        lats = lats[ind]

        
        return(mtopoin, latind, lonind)

    else:
        return(topoin, lats, lons)

def etopo1_subset(file='etopo1.nc', region=None):
    """ read in ardemV2 topography/bathymetry. """
    
    file='/Volumes/WDC_internal/Users/bell/in_and_outbox/MapGrids/etopo_subsets/etopo1_chukchi.nc'
    bathydata = Dataset(file)
    
    topoin = bathydata.variables['Band1'][:]
    lons = bathydata.variables['lon'][:]
    lats = bathydata.variables['lat'][:]
    bathydata.close()
    
    if region == 'Chukchi':
        latbounds = [60,90]
        lonbounds = [-180, -120]

        latind = np.ma.masked_outside(lats,latbounds[0],latbounds[1])
        lonind = np.ma.masked_outside(lons,lonbounds[0],lonbounds[1])
        ind = (latind.mask==True) &(latind.mask == lonind.mask)
        mtopoin = np.ma.masked_array(topoin,ind)
        lons = lons[ind]
        lats = lats[ind]

        
        return(mtopoin, latind, lonind)

    else:
        return(topoin, lats, lons)

"""---------------------------------- Main --------------------------------------------"""

parser = argparse.ArgumentParser(description='Take Lat/Lon pair from excel workbook and find bathydepth from chosen archive ')

parser.add_argument('ExcelDataPath', metavar='ExcelDataPath', type=str, help='Full directory path to excel data')
parser.add_argument('ExcelSheet', metavar='ExcelSheet', type=int, help='Relevant Sheet number in workbook')
parser.add_argument('BathySource', metavar='BathySource', type=str, help='ETOPO5, ARDEMV2, ETOPO1')

args = parser.parse_args()

###Load chosen bathymetry file
if args.BathySource == 'ETOPO5':
	print "Reading Bathy Data from ETOPO5"
	(topoin, lats, lons) = etopo5_data(region='Chukchi')
elif args.BathySource == 'ARDEMV2':
	print "Reading Bathy Data from ARDEM V2"
	(topoin, lats, lons) = ardemV2_data()
elif args.BathySource == 'ETOPO1':
    print "Reading Bathy Data from ETOPO1"
    (topoin, lats, lons) = etopo1_subset()

print "Reading position data from excel"
wb = pd.read_excel(args.ExcelDataPath)

for i, row in wb.iterrows():
    #postive north, negative west
    print "Finding nearest station to {Latitude}, {Longitude}".format(**row)
    model_stn = sphered.nearest_point([row['Latitude'],row['Longitude']],lats,lons, '1d')
    #model_stn = sphered.nearest_point([row['Latitude'],row['Longitude']],lats.data,lons.data, '1d')
    depth = topoin[model_stn[3],model_stn[4]]
    print "{depth}, {topo_lat}, {topo_lon}".format(depth=depth,topo_lat=model_stn[1], topo_lon=model_stn[2])
