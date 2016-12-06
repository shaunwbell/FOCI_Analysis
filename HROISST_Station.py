#!/usr/bin/env

"""
 HROISST_Station.py
 
 Retrieve NCEP HighRes OI SST for defined station:

 Save in EPIC NetCDF standard

 Modified 2016-12-06: Use sbell unified package subroutines for netcdf/time - cleanup code
 
"""
#System Stack
import datetime
import sys
import argparse

#Science Stack
import numpy as np
from netCDF4 import Dataset

# Visual Stack
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# User Stack
from calc import haversine as sphered
from io_utils.EcoFOCI_netCDF_read import EcoFOCI_netCDF
from calc.EPIC2Datetime import EPIC2Datetime, to_UDUNITS, Datetime2EPIC
from utilities import ncutilities as ncutils

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NCEP','Unimak', 'Shumagin','3hr filtered', 'U,V','Winds', 'Gulf of Alaska'

"""------------------------EPIC Write   Modules-------------------------------------------"""


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

def write2epic_cf( file_name, stationid, time, lat_lon, data ):
        ncinstance = ncutil.EPIC_NC_SST_cf(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time=time, latitude=lat_lon[0], longitude=-1 * lat_lon[1], \
            depth_level=0. )
        ncinstance.add_data('T_25', data[0])
        ncinstance.close()
        
"""------------------------- Topo   Modules -------------------------------------------"""

def etopo5_data():
    """ read in etopo5 topography/bathymetry. """
    file = '/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/data/etopo5.nc'
    etopodata = Dataset(file)
    
    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()
    
    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0
    
    lons, lats = np.meshgrid(lons, lats)
    
    return(topoin, lats, lons)
        
"""------------------------- Main   Modules -------------------------------------------"""

parser = argparse.ArgumentParser(description='NCEPHROISST from Single Station')
parser.add_argument('MooringID', metavar='MooringID', type=str, help='MooringID Name')               
parser.add_argument('latitude', metavar='latitude', type=float, help='latitude (+N)')               
parser.add_argument('longitude', metavar='longitude', type=float, help='longitude (+W)')               
parser.add_argument('years', nargs='+', type=int, help='start and stop year')
parser.add_argument('--DataPath', metavar='DataPath', type=str, help='full path to alternate file')
parser.add_argument("-scf",'--store_cf', action="store_true", help='cf conventions - primarily in time')
parser.add_argument("-sep",'--store_epic', action="store_true", help='epic conventions - primarily in time')
parser.add_argument("-plot",'--plot', action="store_true", help='create plot of location')
args = parser.parse_args()

####### CMD Line Parse
### CMD Line User defined Station parameters
station_name = [args.MooringID]
sta_lat = [args.latitude]
sta_long = [args.longitude]

### Hard coded path
if args.DataPath:
    NCEP = args.DataPath
else:
    NCEP = '/Volumes/WDC_internal/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/'

infile = [NCEP + 'sst.day.mean.1981.v2.nc']
print infile

#################
### Grab grid points for future slicing 
#   - assume grid is same in all model output
df = EcoFOCI_netCDF(infile[0])
vars_dic = df.get_vars()
nchandle = df._getnchandle_()
lat_lon = {}
for j, v in enumerate(['lat', 'lon']):
    lat_lon[v] = nchandle.variables[v][:]
df.close()


### Find model points
#Find NCEP nearest point to moorings - haversine formula
# NCEP data is 0->360 (positive east), Moorings are usually expressed +W for FOCI
stn1_pt = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '1d')
stn1_modelpt = [lat_lon['lat'][stn1_pt[3]],lat_lon['lon'][stn1_pt[4]]]

print "stn1 nearest point to {sta_lat}, {sta_lon} which is lat:{sta_modellat} , lon:{sta_modellon}".format(
    sta_lat=sta_lat[0], sta_lon=sta_long[0], 
    sta_modellat=stn1_modelpt[0], sta_modellon=stn1_modelpt[1])

stn1_modelpt[1] = -1.*((180 - stn1_modelpt[1]) + 180)
print "thus converting lon to degrees W positive {sta_modellon}".format(sta_modellon=stn1_modelpt[1])

#loop over all requested data   
years = range(args.years[0],args.years[1]+1)

for yy in years:
    # retrieve only these location's data

    ### SST files
    infile = NCEP + 'sst.day.mean.'+ str(yy) + '.v2.nc'
    print "Working on file " + infile
    df = EcoFOCI_netCDF(infile)
    vars_dic = df.get_vars()
    nchandle = df._getnchandle_()
    print "Parameters availabile: {params}".format(params=vars_dic.keys())
    stn1_data = {}
    for j, v in enumerate(vars_dic): 
        try: #check for nc variable
                stn1_data[v] = nchandle.variables[v][:,stn1_pt[3], stn1_pt[4]]
        except ValueError: #if parameter is not of expected dimensions
            stn1_data[v] = nchandle.variables[v][:]
    stn1_sst = stn1_data['sst']
    df.close()


    #convert to EPIC time
    pydate = to_UDUNITS(stn1_data['time'], 
            time_since_str='days since 1800-01-01 00:00:00')
    epic_time, epic_time1 = Datetime2EPIC(pydat.tolist())

    
    if args.store_epic:
        # write to NetCDF
        outfile = 'data/NOAA_OI_SST_V2_' + args.MooringID + '_' + str(yy) + '.nc'
        print "Writing to Epic NetCDF " + outfile
        write2epic( outfile, station_name[0], [epic_time, epic_time1], stn1_modelpt, [stn1_sst,])
        if args.cf:    
            #days since 1800-1-1 00:00:0.0
            date_str_cf = []
            write2epic_cf( outfile, station_name[0], date_str_cf, stn1_modelpt, [stn1_sst,])
    if args.store_cf:
        # write to NetCDF
        outfile = 'data/NOAA_OI_SST_V2_' + args.MooringID + '_' + str(yy) + '_cf.nc'
        print "Writing to Epic NetCDF " + outfile
        write2epic_cf( outfile, station_name[0], date_str_cf, stn1_modelpt, [stn1_sst,])
        if args.cf:    
            #days since 1800-1-1 00:00:0.0
            date_str_cf = []
            write2epic_cf( outfile, station_name[0], date_str_cf, stn1_modelpt, [stn1_sst,])


if args.plot:
    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=52, \
        urcrnrlat=58,llcrnrlon=-140,urcrnrlon=-130, lat_ts=45)

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
    m.drawmeridians(np.arange(-140,-130,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]) )

    plt.savefig('images/'+args.MoorindID+'.png', bbox_inches='tight', dpi = (100))
    plt.close()


    


