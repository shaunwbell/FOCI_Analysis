#!/usr/bin/env

"""
NARRuv_quiver_plot.py


"""
#System Stack
import datetime

#Science Stack
from netCDF4 import MFDataset, num2date
import numpy as np

# Visual Stack
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter
import matplotlib.ticker as ticker

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2016, 6, 20)
__modified__ = datetime.datetime(2016, 6, 20)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'quiver plot'

"""-------------------------- Initialization params -----------------------------------------"""

### some mpl specif settings for fonts and plot style
mpl.rcParams['svg.fonttype'] = 'none'
plt.style.use('bmh')

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

"""--------------------------------main Routines---------------------------------------"""

####The following is designed to plot U/V vecotrs from WPAK as a stick plot
### Read CF style netcdf files (one time word) from mulitple files
ncfiles = '/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/bell/jul/'
nchandle = MFDataset(ncfiles+'*thin.nc',aggdim="time")
vars_dic = get_vars(nchandle)
data1 = ncreadfile_dic(nchandle,vars_dic.keys())
nchandle.close()

### some data manipulation or massaging for plots
## generate daily averages from hourly data and label it 12Z of the day averaged
subset_t, ucomp, vcomp = [], [], []
'''
total_date_range = np.arange(data1["time"].min(),data1["time"].max(),1)
for count in total_date_range:
    tind = (np.where(data1["time"] % count < 1))
    ucomp = np.hstack((ucomp, np.mean(data1['WU_422'][tind,0,0,0])))
    vcomp = np.hstack((vcomp, np.mean(data1['WV_423'][tind,0,0,0])))
    subset_t = np.hstack((subset_t, count+0.5))
xdate = num2date(subset_t, "days since 1800-01-01")
xdate = [x.toordinal() for x in xdate]
'''
xdate=[(num2date(x, "days since 1800-01-01")).hour/24. + (num2date(x, "days since 1800-01-01")).toordinal() for x in data1["time"]]
ucomp=data1['WU_422'][:,0,0,0]
vcomp=data1['WV_423'][:,0,0,0]
#exchange 1e35 for 0
ucomp[ucomp >= 1e30] = np.nan
vcomp[vcomp >= 1e30] = np.nan

magnitude = np.sqrt(ucomp**2 + vcomp**2)


####The following is designed to plot U/V vecotrs from NARR as a stick plot
### Read CF style netcdf files (one time word) from mulitple files
ncfiles = '/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/data/'
nchandle = MFDataset(ncfiles+'*cf.nc',aggdim="time")
vars_dic = get_vars(nchandle)
data2 = ncreadfile_dic(nchandle,vars_dic.keys())
nchandle.close()

### some data manipulation or massaging for plots
## generate daily averages from hourly data and label it 12Z of the day averaged
subset_t, ucomp_narr, vcomp_narr = [], [], []
'''
total_date_range = np.arange(data2["time"].min(),data2["time"].max(),1)
for count in total_date_range:
    tind = (np.where(data2["time"] % count < 1))
    ucomp_narr = np.hstack((ucomp_narr, np.mean(data2['WU_422'][tind,0,0,0])))
    vcomp_narr = np.hstack((vcomp_narr, np.mean(data2['WV_423'][tind,0,0,0])))
    subset_t = np.hstack((subset_t, count+0.5))
xdate_narr = num2date(subset_t, "days since 1800-01-01")
xdate_narr = [x.toordinal() for x in xdate_narr]
'''
xdate_narr=[(num2date(x, "days since 1900-01-01")).hour/24. + (num2date(x, "days since 1900-01-01")).toordinal() for x in data2["time"]]
ucomp_narr=data2['WU_422'][:,0,0,0]
vcomp_narr=data2['WV_423'][:,0,0,0]
#exchange 1e35 for 0
ucomp_narr[ucomp_narr >= 1e30] = np.nan
vcomp_narr[vcomp_narr >= 1e30] = np.nan



### Quiver / Stick plot
# Plot quiver
fig = plt.figure(1)
for i in range(1,12):
    ax1 = fig.add_subplot(11,1,i)
    # 1D Quiver plot
 
    q = ax1.quiver(xdate,0,ucomp,vcomp,color='r',units='y',scale_units='y',
                   scale = 1,headlength=2,headaxislength=2,width=0.1,alpha=.95)
    qk = plt.quiverkey(q,0.2, 0.05, 5,r'$5 \frac{m}{s}$',labelpos='W',
                   fontproperties={'weight': 'bold'})
    q = ax1.quiver(xdate_narr,0,ucomp_narr,vcomp_narr,color='k',units='y',scale_units='y',
                   scale = 1,headlength=2,headaxislength=2,width=0.1,alpha=.95)
    qk = plt.quiverkey(q,0.2, 0.05, 5,r'$5 \frac{m}{s}$',labelpos='W',
                   fontproperties={'weight': 'bold'}) 
    ax1.set_ylim(np.nanmin(vcomp), np.nanmax(vcomp))
    ax1.set_ylabel("(m/s)")
    ax1.set_xlim(datetime.datetime(2005+i-1, 01, 01).toordinal(),datetime.datetime(2005+i-1, 12, 31).toordinal())
    ax1.xaxis.set_major_locator(MonthLocator())
    ax1.xaxis.set_minor_locator(MonthLocator(bymonth=range(1,13), bymonthday=15))
    ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    ax1.xaxis.set_minor_formatter(DateFormatter('%b %y'))
    ax1.tick_params(axis='both', which='minor', labelsize=10)

t = fig.suptitle('WPAK Daily Averaged Wind', fontsize=12)
t.set_y(0.03)
fig.autofmt_xdate()
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*6) )
plt.savefig('images/M2_wpakwinds.png', bbox_inches='tight', dpi = (300))
plt.close()
#plt.show()
