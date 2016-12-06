#!/usr/bin/env

"""
 GOA_Winds_StormPatterns.py
 
 Compare Gorepoint/globec winds (along shore) to telleconnection indices
    GorePoint - 58deg 58min N, 150deg 56min W 
    and Globec3 59.273701N, 148.9653W

 Files are created by GOA_Winds_NARR_model_prep.py
 
 -Filtered NARR winds with a triangular filter (1/4, 1/2, 1/4) and output every 3hrs
 -Provided U, V 

 -Saved in EPIC NetCDF standard
"""
#System Stack
import datetime
import sys

#Science Stack
import numpy as np
from scipy import stats

# User Stack
from utilities import ncutilities as ncutil

# Visual Stack
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 04, 29)
__modified__ = datetime.datetime(2014, 04, 29)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'NARR','GLOBEC3', 'Gorept','AO/NAO/PNA', 'U,V','Winds', 'Gulf of Alaska'

"""------------------------General   Modules-------------------------------------------"""

def from_netcdf(infile):
    """ Uses ncreadfile_dic which returns a dictionary of all data from netcdf"""

    ###nc readin/out
    nchandle = ncutil.ncopen(infile)
    params = ncutil.get_vars(nchandle) #gets all of them


    ncdata = ncutil.ncreadfile_dic(nchandle, params)
    ncutil.ncclose(nchandle)
    
    return (ncdata, params)
    
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
    elif file_flag == 'Index':
        """ yyyy mm dd"""
        python_time=datetime.datetime.strptime(file_time,'%Y %m %d').toordinal()
    else:
        print "time flag not recognized"
        sys.exit()
        
    return python_time


"""------------------------- MATH   Modules -------------------------------------------"""
def hourly_2_ave(ltbound,utbound, time, data, time_base=6.):
    """ bin average times into specified bins """
    interval = time_base / 24.
    tarray = np.arange(ltbound, utbound,interval)
    dmean = np.zeros_like(tarray) * np.nan
    dstd = np.zeros_like(tarray) * np.nan

    for i, val in enumerate(tarray):
        ind = (time >= val) & (time < val+interval)
        dmean[i] = data[ind].mean()
        dstd[i] = data[ind].std()
    
    return  { 'dtime':tarray, 'dmean':dmean ,'dstd':dstd,} 
    
def rotate_coord(angle_rot, mag, direct):
    """ converts math coords to along/cross shelf.
    + onshore  / along coast with land to right (right handed)
    - offshore / along coast with land to left
    
    Todo: convert met standard for winds (left handed coordinate system
    """
    
    direct =  direct - angle_rot
    
    along = mag * np.sin(np.deg2rad(direct))
    cross = mag * np.cos(np.deg2rad(direct))
    
    return (along, cross)    

def lin_fit(x, y):
    """ scipy linear regression routine"""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

    return ( slope, intercept, r_value, p_value, std_err )

def moving_average(x, n, type='simple'):
    """
    compute an n period moving average.

    type is 'simple' | 'exponential'

    """
    #x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()


    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a
        
"""------------------------- Main   Modules -------------------------------------------"""

### READ AO/PNA indices from txt file
AO_file = '/Users/bell/Data_Local/teleconnections/norm.daily.ao.index.b500101.current.ascii'
PNA_file = '/Users/bell/Data_Local/teleconnections/norm.daily.pna.index.b500101.current.ascii'

# ingest indicies
PNA_index, PNA_time = [], [] #some missing ind
with open(PNA_file, 'rb') as f:
    for k, line in enumerate(f.readlines()):
        PNA_index = PNA_index + [line.strip().split()[-1]]
        PNA_time = PNA_time + [date2pydate(" ".join(line.strip().split()[:-1]), file_flag='Index')]
        
PNA_index= np.array(PNA_index, float)
PNA_time= np.array(PNA_time)

AO_index, AO_time = [], [] #some missing ind
with open(AO_file, 'rb') as f:
    for k, line in enumerate(f.readlines()):
        AO_index = AO_index + [line.strip().split()[-1]]
        AO_time = AO_time + [date2pydate(" ".join(line.strip().split()[:-1]), file_flag='Index')]
        
AO_index= np.array(AO_index, float)
AO_time= np.array(AO_time)

### NARR wind files (preprocessed) for specific locations - winds have a triangle filter on them
NARR = '/Users/bell/Programs/Python/FOCI_Analysis/GOA_Winds/data/'

station_name = ['Globec3','GorePt']
sta_lat = [59.273701,58.9666666666666667]
sta_long = [148.9653,150.9333333333333333]

#loop over all requested data    
years = range(1984, 2014,1)

NARR_time = []
NARR_uwnd = []
NARR_vwnd = []
for iyear in years:
    globec3_data, NARRkeys = from_netcdf(NARR+'NARR_globec_'+str(iyear)+'.nc')
    NARR_time = NARR_time + date2pydate(globec3_data['time'], globec3_data['time2'])
    NARR_uwnd = np.append(NARR_uwnd, globec3_data['WU_422'][:,0,0,0])
    NARR_vwnd = np.append(NARR_vwnd, globec3_data['WV_423'][:,0,0,0])    

NARR_time = np.array(NARR_time)

### daily averages
time_bin = 24.

NARRDaily_uwnd = hourly_2_ave(NARR_time.min(),NARR_time.max(), NARR_time, NARR_uwnd, time_base=time_bin)    
NARRDaily_vwnd = hourly_2_ave(NARR_time.min(),NARR_time.max(), NARR_time, NARR_vwnd, time_base=time_bin)    

NARR_wndmag = np.sqrt((NARRDaily_uwnd['dmean']**2)+(NARRDaily_vwnd['dmean']**2))
NARR_wind_dir_math = np.rad2deg(np.arctan2(NARRDaily_vwnd['dmean'] , NARRDaily_uwnd['dmean']))
NARR_along,NARR_across = rotate_coord(120., NARR_wndmag, NARR_wind_dir_math)

"""----"""
# Calculate correlations for 3month spans
corr_PNA = {}
corr_AO = {}
for drange in range(1980,2014,1):
    for mrange in range(1,12,3):
        start_ind = datetime.datetime.strptime(str(drange) + ' ' + str(mrange) + ' 01','%Y %m %d').toordinal()
        end_ind = datetime.datetime.strptime(str(drange) + ' ' + str(mrange+2) + ' 01','%Y %m %d').toordinal()
        PNA_ind = (PNA_time >= start_ind) & (PNA_time <= end_ind)
        AO_ind = (AO_time >= start_ind) & (AO_time <= end_ind)
        NARR_ind = (NARRDaily_uwnd['dtime'] >= start_ind) & (NARRDaily_uwnd['dtime'] <= end_ind)
#        NARR_along_stand = (NARR_along - np.nanmean(NARR_along)) / np.nanstd(NARR_along)
#        NARR_along_stand = (NARR_along - np.nanmin(NARR_along)) / (np.nanmax(NARR_along) - np.nanmin(NARR_along)) #actually normalized - ignor var name
#        PNA_index_stand = (PNA_index - np.nanmean(PNA_index)) / (np.nanstd(PNA_index)) 
#        PNA_index_stand = (PNA_index - np.nanmin(PNA_index)) / (np.nanmax(PNA_index) - np.nanmin(PNA_index)) #actually normalized - ignor var name
        if not np.size(NARR_along[NARR_ind]) == 0:
            #(slope, intercept, r_value, p_value, std_err) = lin_fit(NARR_along_stand[NARR_ind], PNA_index[PNA_ind])
            #corr[start_ind] = r_value**2
            corr_PNA[start_ind] = np.corrcoef(NARR_along[NARR_ind], PNA_index[PNA_ind])[0][1]
            corr_AO[start_ind] = np.corrcoef(NARR_along[NARR_ind], AO_index[PNA_ind])[0][1]
        else:
            corr_PNA[start_ind] = 0.0
            corr_AO[start_ind] = 0.0

# 30day running filter for wind timeseries 
NARR_along_rm = moving_average(NARR_along,(30))
AO_index_rm = moving_average(AO_index,(30))
PNA_index_rm = moving_average(PNA_index,(30))
"""------------------------- Plotting Modules -------------------------------------------"""

year_bounds = [[datetime.datetime.strptime('1980 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('1985 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('1990 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('1995 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2000 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2005 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2010 01 01','%Y %m %d').toordinal()],
                [datetime.datetime.strptime('1985 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('1990 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('1995 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2000 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2005 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2010 01 01','%Y %m %d').toordinal(),
                    datetime.datetime.strptime('2015 01 01','%Y %m %d').toordinal()]]
fig = plt.figure()
for splot in range(0,7,1):
    ax1 = plt.subplot(7,1,splot+1)
    plt.plot(NARRDaily_uwnd['dtime'], NARR_along_rm, 'r')
    for i,kk in enumerate(corr_PNA.keys()):
        
        if (corr_PNA[kk]) >=.2 and (corr_PNA[kk]) <=.3:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.1)
        elif (corr_PNA[kk]) >=.3 and (corr_PNA[kk]) <=.4:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.3)
        elif (corr_PNA[kk]) >=.4 and (corr_PNA[kk]) <=.5:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.5)
        elif (corr_PNA[kk]) >=.5:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.8)
    
        elif (corr_PNA[kk]) <=-0.2 and (corr_PNA[kk]) >=-0.3:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.1)
        elif (corr_PNA[kk]) <=-0.3 and (corr_PNA[kk]) >=-0.4:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.3)
        elif (corr_PNA[kk]) <=-0.4 and (corr_PNA[kk]) >=-0.5:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.5)
        elif (corr_PNA[kk]) <=-0.5:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.8)
            
    ax1.set_ylim((-10,10))
    ax2 = ax1.twinx()
    plt.plot(PNA_time, PNA_index_rm, 'b')
    ax2.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax2.set_ylim((-3,3))
    ax2.set_xlim(year_bounds[0][splot],year_bounds[1][splot])
    ax2.xaxis.set_major_locator(MonthLocator(bymonth=[3,10], bymonthday=1))

fig.suptitle('NARR Along-Shore Winds corr PNA Index at Globec3')

DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )
plt.savefig('NARR_along_PNA_globec.png', bbox_inches='tight', dpi = (100))
plt.close()
    
    
fig = plt.figure()
for splot in range(0,7,1):
    ax1 = plt.subplot(7,1,splot+1)
    plt.plot(NARRDaily_uwnd['dtime'], NARR_along_rm, 'r')
    for i,kk in enumerate(corr_AO.keys()):
        
        if (corr_AO[kk]) >=.2 and (corr_AO[kk]) <=.3:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.1)
        elif (corr_AO[kk]) >=.3 and (corr_AO[kk]) <=.4:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.3)
        elif (corr_AO[kk]) >=.4 and (corr_AO[kk]) <=.5:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.5)
        elif (corr_AO[kk]) >=.5:
            fill1 = ax1.axvspan(kk, kk+90, color='k', alpha=0.8)
    
        elif (corr_AO[kk]) <=-0.2 and (corr_AO[kk]) >=-0.3:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.1)
        elif (corr_AO[kk]) <=-0.3 and (corr_AO[kk]) >=-0.4:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.3)
        elif (corr_AO[kk]) <=-0.4 and (corr_AO[kk]) >=-0.5:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.5)
        elif (corr_AO[kk]) <=-0.5:
            fill1 = ax1.axvspan(kk, kk+90, color='r', alpha=0.8)
            
    ax1.set_ylim((-10,10))
    ax2 = ax1.twinx()
    plt.plot(AO_time, AO_index_rm, 'b')
    ax2.xaxis.set_major_formatter(DateFormatter('%b %Y'))
    ax2.set_ylim((-3,3))
    ax2.set_xlim(year_bounds[0][splot],year_bounds[1][splot])
    ax2.xaxis.set_major_locator(MonthLocator(bymonth=[3,10], bymonthday=1))

fig.suptitle('NARR Along-Shore Winds corr AO Index at Globec3')
    
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )
plt.savefig('NARR_along_AO_globec.png', bbox_inches='tight', dpi = (100))
plt.close()