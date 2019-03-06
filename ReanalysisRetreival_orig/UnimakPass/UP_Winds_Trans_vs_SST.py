#!/usr/bin/env

"""
 UP_Winds_Trans_vs_SST.py
 
    Using U,V (6hr from NARR) to calculate a transport index
    
    Using SST (daily) from HR
    
 NARR U/V winds (triangel filtered and subsampled to 6 hours)
 ----
 NCEP Reanalysis data provided by the NOAA/OAR/ESRL PSD, Boulder,
  Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/ 
 
 SST -
 ---
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
import pandas as pd

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
"""--------------------------------time Routines---------------------------------------"""

def date2pydate(file_time, file_time2=None, time_since_str=None, file_flag='EPIC'):


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
    elif file_flag == 'netCDF4':
        """ Use time_since_str"""
        python_time = num2date(file_time,'time_since_str','standard')
    else:
        print "time flag not recognized"
        sys.exit()
        
    return np.array(python_time)

"""--------------------------------netcdf Routines---------------------------------------"""
def get_global_atts(nchandle):

    g_atts = {}
    att_names = nchandle.ncattrs()
    
    for name in att_names:
        g_atts[name] = nchandle.getncattr(name)
        
    return g_atts

def get_vars(nchandle):
    return nchandle.variables

def ncreadfile_dic(nchandle, params):
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
                data[v] = nchandle.variables[v][:]

        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)
   
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
    

"""------------------------- Main   Modules -------------------------------------------"""

### list of files
NARR_dir = '/Users/bell/in_and_outbox/2016/stabeno/feb/unimakwinds_narr/shumigan_downstream/'
HROISST_dir = '/Users/bell/in_and_outbox/2016/stabeno/feb/unimakwinds_narr/shumigan_downstream_sst/'

#loop over every year from 1981 to 2015.
# calculate desired average (based on time stamps)
data_flag = 'winds'
for year in range(2016,2019):
    #print "Working on year {0}".format(year)
    sstfile = HROISST_dir + 'NOAA_OI_SST_V2_stn1_' + str(year) + '.nc'
    uvfile = NARR_dir + 'NARR_stn1_' + str(year) + '.nc'
    
    if data_flag == 'sst':
        #open sst file
        ###nc readin
        nchandle = Dataset(sstfile,'a')
        global_atts = get_global_atts(nchandle)
        vars_dic = get_vars(nchandle)
        sstdata = ncreadfile_dic(nchandle, vars_dic.keys())
        nchandle.close()
    
        ## Create and save monthly mean data of SST
        tmptime = date2pydate(sstdata['time'],sstdata['time2'])
        ssttime_month = [datetime.datetime.fromordinal(int(x)).month for x in tmptime]
        for month in range(1,13):
            tind = np.where(np.array(ssttime_month) == month)
            sstmean = np.mean(sstdata['T_25'][tind,0,0,0])
            print "{0}-{1}-01, {2}".format(year,month,sstmean)            
            
    elif data_flag == 'winds':        
        #open uv file
        ###nc readin
        nchandle = Dataset(uvfile,'a')
        global_atts = get_global_atts(nchandle)
        vars_dic = get_vars(nchandle)
        uvdata = ncreadfile_dic(nchandle, vars_dic.keys())
        nchandle.close()
    
        ## Create and save monthly mean data of UV winds/transport
        tmptime = date2pydate(uvdata['time'],uvdata['time2'])
        uvtime_month = [datetime.datetime.fromordinal(int(x)).month for x in tmptime]
        for month in range(1,13):
            tind = np.where(np.array(uvtime_month) == month)
            umean = np.mean(uvdata['WU_422'][tind,0,0,0])    
            vmean = np.mean(uvdata['WV_423'][tind,0,0,0])    
            print "{0}-{1}-01, {2}, {3}".format(year,month,umean,vmean)            
    else:
        print "skipping data read"

plot_flag = False
#After all data has been averaged plot
#plot as timeseries but colorcode the values as follows:
#  Break sst data into quintiles (first hack, find max and min and evenly divide by 5)
#  Colorcode transport as a function of sst quintiles where warm 1/5 is bright red, cold 1/5 is bright blue
#  2/5's are much lighter and median/mean is grey
if plot_flag:
    #watch for nan's with extra spaces
    data = pd.read_csv('UP_transport.csv')
    maxd = data.max()['SST']
    mind = data.min()['SST']
    bounds = np.arange(mind,maxd,(maxd-mind)/5)
    bounds = np.hstack([bounds,maxd])
    
    mag =np.sqrt(data['U']**2 + data['V']**2)
    ang = np.rad2deg(np.arctan2(data['V'],data['U']))
    transport_rough = rotate_coord(45, mag, ang)
    print transport_rough[0]
    transport_rough = mag*np.sqrt(2)/2
    for ind,val in enumerate(bounds):
        if ind == 0:
            pind = np.where((data['SST'] >= val) & (data['SST'] <= bounds[ind+1]))
            plt.plot(pd.to_datetime(data['Date'][pind[0]]),transport_rough[pind[0]],'.',color='#001AFF',markersize=20)
        if ind == 1:
            pind = np.where((data['SST'] >= val) & (data['SST'] <= bounds[ind+1]))
            plt.plot(pd.to_datetime(data['Date'][pind[0]]),transport_rough[pind[0]],'.',color='#7D8AFF',markersize=20)
        if ind == 2:
            pind = np.where((data['SST'] >= val) & (data['SST'] <= bounds[ind+1]))
            plt.plot(pd.to_datetime(data['Date'][pind[0]]),transport_rough[pind[0]],'.',color='#B2B2B2',markersize=20)
        if ind == 3:
            pind = np.where((data['SST'] >= val) & (data['SST'] <= bounds[ind+1]))
            plt.plot(pd.to_datetime(data['Date'][pind[0]]),transport_rough[pind[0]],'.',color='#FF9B9B',markersize=20)
        if ind == 4:
            pind = np.where((data['SST'] >= val) & (data['SST'] <= bounds[ind+1]))
            plt.plot(pd.to_datetime(data['Date'][pind[0]]),transport_rough[pind[0]],'.',color='#FF0000',markersize=20)
        if ind == 5:
            break    