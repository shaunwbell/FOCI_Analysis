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
import argparse

#Science Stack
import numpy as np
from netCDF4 import Dataset

# User Stack
import utilities.haversine as sphered
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
    date, time, uwnd, vwnd, atemp, bpress = [], [], [], [], [], []
    with open(ifile, 'rb') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader) #skip header
        """   DAT   TIM          WU           WV           AT           BP        """
        for row in csv_reader:
            try:
                r0,r1,r2,r3,r4,r5,r6 = row[0].strip().split()
            except ValueError:
                r0,r1,r2,r3,r4,r5 = row[0].strip().split()
            date.append(r0)
            time.append(r1)
            uwnd.append(r2)
            vwnd.append(r3)


    return {'DAT': np.array(date, int), 'TIM':np.array(time, float), 'WU':np.array(uwnd, float),\
     'WV':np.array(vwnd, float)}
     
def write2epic( file_name, stationid, time, lat_lon, data ):
        ncinstance = ncutil.EPIC_NC_rad(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time[0]))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time1=time[0], time2=time[1], latitude=lat_lon[0], longitude=-1 * lat_lon[1], \
            depth_level=10. )
        ncinstance.add_data('Qs_133', data[0])

        ncinstance.close()

def write2epic_cf( file_name, stationid, time, lat_lon, data ):
        ncinstance = ncutil.EPIC_NC_SST_cf(savefile=file_name)
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts(Station_Name=stationid, file_name=( __file__.split('/')[-1]) )
        ncinstance.dimension_init(len_time=len(time))
        ncinstance.variable_init()
        ncinstance.add_coord_data(time=time, latitude=lat_lon[0], longitude=-1 * lat_lon[1], \
            depth_level=10. )
        ncinstance.add_data('Qs_133', data[0])

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

"---"

def triangle_smoothing(data_in):
    weights=np.array([0.25,0.5,0.25])
    filtered_data = np.convolve(data_in,np.array(weights),'same') #edge effects
    
    return filtered_data

"""------------------------- Topo   Modules -------------------------------------------"""

def etopo5_data():
    """ read in etopo5 topography/bathymetry. """
    file = '/Volumes/WDC_internal/Users/bell/in_and_outbox/MapGrids/etopo5.nc'
    etopodata = Dataset(file)
    
    topoin = etopodata.variables['bath'][:]
    lons = etopodata.variables['X'][:]
    lats = etopodata.variables['Y'][:]
    etopodata.close()
    
    topoin,lons = shiftgrid(0.,topoin,lons,start=False) # -360 -> 0
    
    lons, lats = np.meshgrid(lons, lats)
    
    return(topoin, lats, lons)
        
"""------------------------- Main   Modules -------------------------------------------"""

parser = argparse.ArgumentParser(description='NARR from Single Station')
parser.add_argument('MooringID', metavar='MooringID', type=str, help='MooringID Name')               
parser.add_argument('latitude', metavar='latitude', type=float, help='latitude (+N)')               
parser.add_argument('longitude', metavar='longitude', type=float, help='longitude (+W)')               
parser.add_argument('years', nargs='+', type=int, help='start and stop year')
parser.add_argument('--DataPath', metavar='DataPath', type=str, help='full path to alternate file')
parser.add_argument("-cf",'--cf', action="store_true", help='cf conventions - primarily in time')
parser.add_argument("-sm",'--sitemap', action="store_true", help='create a site map')
args = parser.parse_args()


### list of files
if args.DataPath:
    NARR = args.DataPath
else:
    NARR = '/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/'

ftype = 'lhtfl'
infile = [NARR + ftype + '.2017.nc'] #used just to get grid sections

print infile
### Grab grid points for future slicing - assume grid is same in all model output
lat_lon = latlon_grid(infile[0])

station_name = [args.MooringID]
sta_lat = [args.latitude]
sta_long = [args.longitude]

#Find NARR nearest point to moorings - haversine formula
station_1 = sphered.nearest_point([sta_lat[0],-1 * sta_long[0]],lat_lon['lat'],lat_lon['lon'], '2d')
station_1_modelpt = [lat_lon['lat'][station_1[3],station_1[4]],lat_lon['lon'][station_1[3],station_1[4]]]

print "station_1 nearest point to %s, %s which is lat:%s , lon:%s" \
    % (sta_lat[0], sta_long[0], station_1_modelpt[0], station_1_modelpt[1])
    
#loop over all requested data   
years = range(args.years[0],args.years[1]+1)

for yy in years:
    # retrieve only these location's data
    # uwnd
    infile = NARR + ftype+ '.'+ str(yy) + '.nc'
    print "Working on file " + infile
    station_1_data = from_netcdf_1dsplice(infile, None, station_1[3], station_1[4])

    #filter data
    station_1dswrf = station_1_data[ftype]
    
    #convert to EPIC time
    pydate = date2pydate(station_1_data['time'], file_flag='NARR')
    epic_time, epic_time1 = pydate2EPIC(pydate)

    # output u,v wind components from model grid points
    save_to_nc = True
    if save_to_nc:
        # write to NetCDF
        outfile = 'data/NARR_' + station_name[0] + '_' + str(yy) + '_'+ ftype+'.nc'
        print "Writing to Epic NetCDF " + outfile
        if args.cf:    
            #days since 1800-1-1 00:00:0.0
            date_str_cf = []
            write2epic_cf( outfile, station_name[0], date_str_cf, stn1_modelpt, [station_1dswrf])
        else:
            write2epic( outfile, station_name[0], [epic_time, epic_time1], station_1_modelpt, [station_1dswrf])

if args.sitemap:
    (topoin, elats, elons) = etopo5_data()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    m = Basemap(resolution='i',projection='merc', llcrnrlat=55, \
        urcrnrlat=75,llcrnrlon=-180,urcrnrlon=-145, lat_ts=60)
    


    # Mooring Data
    x_moor, y_moor = m([-1. * sta_long[0],],sta_lat)
    x_close, y_close = m([station_1_modelpt[1],], [station_1_modelpt[0],])

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
    m.drawparallels(np.arange(60,75,2.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
    m.drawmeridians(np.arange(-180,-145,2.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
    #m.fillcontinents(color='black')

    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0], DefaultSize[1]) )

    plt.savefig('images/Chuckchi_region.png', bbox_inches='tight', dpi = (100))
    plt.close()
