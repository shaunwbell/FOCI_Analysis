#!/usr/bin/env

"""
Program: SST_Region.py

Purpose: generate a timeseries of SST data at/near EcoFOCI M2 (Berring Sea) Mooring

DataSource: ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/

Input/Arguments: Latitude and Longitude, raw data file
 
"""
#System Stack
import datetime

#Science Stack
from netCDF4 import Dataset
import numpy as np

# User Stack
import utilities.ncutilities as ncutil

#Visual Packages
import matplotlib as mpl
mpl.use('Agg') 
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
            pytime = []
    
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

ncfile_all = ['/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/sst.day.mean.2011.v2.nc', \
          '/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/sst.day.mean.2012.v2.nc', \
          '/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/sst.day.mean.2013.v2.nc', \
          '/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/sst.day.mean.2014.v2.nc', \
          '/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/sst.day.mean.2015.v2.nc', ]

for ncfile in ncfile_all:
    print "Working with {0}".format(ncfile)
    nchandle = Dataset(ncfile,'r')
    params = get_vars(nchandle)
    for ind in range(0,365):
        MapData = ncreadfile_dic_oneday(nchandle,params,ind,['lat','lon'])
        time = date2pydate(MapData['time'], file_flag='NCEP_days')
        time_label = datetime.datetime.fromordinal(time).strftime("%Y-%m-%d")


        """------------------------------------------------------------------------------------"""

        print "Generating image"
        ## plot
        #(topoin_tot, elats_tot, elons_tot) = etopo5_data()
        #(topoin, elats, elons) = etopo5_data()
        (topoin, elats, elons) = (MapData['sst'], MapData['lat'], MapData['lon'])
        topoin,elons = shiftgrid(180.,topoin,elons,start=False)

        #build regional subset of data
        data_bounds = {'latmin': 60, 'latmax':72, 'lonmin': -150, 'lonmax':-180}
        plot_bounds = {'latmin': 62, 'latmax':70, 'lonmin': -155, 'lonmax':-180}

        topoin = topoin[find_nearest(elats,data_bounds['latmin']):find_nearest(elats,data_bounds['latmax']),\
                        find_nearest(elons,data_bounds['lonmax']):find_nearest(elons,data_bounds['lonmin'])]
        elons = elons[find_nearest(elons,data_bounds['lonmax']):find_nearest(elons,data_bounds['lonmin'])]
        elats = elats[find_nearest(elats,data_bounds['latmin']):find_nearest(elats,data_bounds['latmax'])]


        fig = plt.figure()
        ax = plt.subplot(111)
        m = Basemap(resolution='i',projection='merc', llcrnrlat=plot_bounds['latmin'], \
            urcrnrlat=plot_bounds['latmax'],llcrnrlon=plot_bounds['lonmax'],urcrnrlon=plot_bounds['lonmin'],\
            lat_ts=45)


        # Site Data
        #x_moor, y_moor = m(-164.051,56.8678)

        #ETOPO 5 contour data 
        elons, elats = np.meshgrid(elons, elats)
        ex, ey = m(elons, elats)

        #CS = m.imshow(topoin, cmap='Greys_r') #
        CS_l = m.contour(ex,ey,topoin, levels=[0,2,4,6,8,10,12,14,16], linestyle='--', linewidths=1, colors='black', alpha=.75) 
        im1 = m.pcolormesh(ex,ey,topoin,vmin=0,vmax=16,cmap=plt.cm.RdBu_r)
        #CS = m.contourf(ex,ey,topoin, levels=[-8,7,6,5,4,3,2,1,0,]) 
        #plt.clabel(CS_l, inline=1, fontsize=8, fmt='%1.0f')
        cb = m.colorbar(im1,"bottom", size="5%", pad="10%")

        #plot points
        #m.scatter(x_moor,y_moor,16,marker='o',facecolors='none', edgecolors='k', alpha=.25)
        #plt.text(x_moor-10000,y_moor-10000,'M2', fontsize=16 ) 
    
        #m.drawcountries(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(np.arange(46,80,4.),labels=[1,0,0,0],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw parallels
        m.drawmeridians(np.arange(-180,-140,5.),labels=[0,0,0,1],color='black',dashes=[1,1],labelstyle='+/-',linewidth=0.2) # draw meridians
        m.fillcontinents(color='0.3')

        plt.annotate(time_label, xy=(0.025, .025), color='k', xycoords='axes fraction', fontsize=16)

        DefaultSize = fig.get_size_inches()
        fig.set_size_inches( (DefaultSize[0]*1.5, DefaultSize[1]*2) )

        plt.savefig('images/sst/' + time_label + '_sst_map.png', bbox_inches='tight', dpi = (300))
        plt.close()

    nchandle.close()
