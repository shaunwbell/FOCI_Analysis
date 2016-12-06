#!/usr/bin/env

"""
 Chuckchi_Winds_NARR_model_prep.py
 
 Retrieve NARR winds for one locations
 Icy Cape Line, Ckip2
 Latitude = 70.8401 Longitude = 163.2054

 Filter NARR winds with a triangular filter (1/4, 1/2, 1/4) and output every 3hrs
 Provide U, V 

 Save in EPIC NetCDF standard
"""
#System Stack
import datetime
import csv

#Science Stack
import numpy as np

# User Stack
from utilities import ncutilities as ncutil


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 01, 13)
__modified__ = datetime.datetime(2014, 01, 13)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR','station_1','3hr filtered', 'U,V','Winds', 'Chuckchi'

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

def csvread(ifile):
    date, time, uwnd, vwnd = [], [], [], []
    with open(ifile, 'rb') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader) #skip header
        next(csv_reader) #skip header
        next(csv_reader) #skip header
        next(csv_reader) #skip header
        next(csv_reader) #skip header
        next(csv_reader) #skip header
        for row in csv_reader:
            r0,r1,r2,r3,r4,r5 = row[0].strip().split()

            date.append(r0)
            time.append(r1)
            uwnd.append(r4)
            vwnd.append(r5)


    return date, time, uwnd, vwnd
     
def write2epic( file_name, stationid, time, lat_lon, data ):
        ncinstance = ncutil.EPIC_NC(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time[0]))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time1=time[0], time2=time[1], latitude=lat_lon[0], longitude=-1 * lat_lon[1][0], \
            depth_level=10. )
        ncinstance.add_data('WU_422', data[0])
        ncinstance.add_data('WV_423', data[1])
        ncinstance.add_data('AT_21', data[2])
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

        
"""------------------------- Main   Modules -------------------------------------------"""

### list of files
NARR = '/Users/bell/Data_Local/from_phyllis/'
infile = NARR + 'narr_c2.dat' #used just to get grid sections

#(sta1,)
station_name = [ 'C2',]
sta_lat = [71.2,]
sta_long = [164.3,]

#read in file
#01-JAN-2010 00:00 /     1:  -2.10  -2.20
header = 6
date_c2,time_c2,uwnd_c2,vwnd_c2 = csvread(infile)



#convert to EPIC time
pydate , station_u, station_v, station_air = [], [], [], []
for i,v in enumerate(date_c2):
    tempdate = datetime.datetime.strptime(date_c2[i]+' '+time_c2[i],'%d-%b-%Y %H:%M')
    #find only 0,6,12,18 time stamps
    if (tempdate.hour == 0) or (tempdate.hour == 6) or (tempdate.hour == 12) or (tempdate.hour == 18):
        pydate.append(tempdate.toordinal() + tempdate.hour /24.)
        station_u.append(uwnd_c2[i])
        station_v.append(vwnd_c2[i])
        station_air.append(1e35)
        
epic_time, epic_time1 = pydate2EPIC(pydate)

# output u,v wind components from model grid points
save_to_nc = True
if save_to_nc:
    # write to NetCDF
    outfile = 'data/narr_c2.nc'
    print "Writing to Epic NetCDF " + outfile
    write2epic( outfile, station_name[0], [epic_time, epic_time1], [sta_lat, sta_long], [station_u, station_v, station_air])
