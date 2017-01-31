# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:42:45 2017

@author: bell
"""


import glob
import argparse
import datetime

from ftplib import FTP

import xarray as xa

#Visual Packages
import matplotlib as mpl
#mpl.use('Agg') 
import matplotlib.pyplot as plt

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2017, 1, 30)
__modified__ = datetime.datetime(2017, 1, 30)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'xarray','SST'

"""--------------------------------helper Routines---------------------------------------"""


"""--------------------------------main Routines---------------------------------------"""

parser = argparse.ArgumentParser(description='NASA WorldView MODIS Retrieval')
parser.add_argument('year', metavar='year', type=str, help='year of data set')               
parser.add_argument('dataset', metavar='dataset', type=str, help='dataset ("anom","err","mean")' )
parser.add_argument('-plot','--plot', action="store_true", help='plot last 4 weeks' )
                  
args = parser.parse_args()


ftp = FTP('ftp.cdc.noaa.gov',user='anonymous')
ftp.cwd('Datasets/noaa.oisst.v2.highres')
path = '/Volumes/WDC_internal/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/'
filename = 'sst.day.'+args.dataset+'.'+args.year+'.v2.nc'
# download
fhandle = open(path + filename, 'wb')
ftp.retrbinary('RETR ' + filename, fhandle.write)                                       # Imprimimos por pantalla lo que estamos descargando        #fhandle.close()
fhandle.close()  
ftp.close()

if args.plot:
	df = xa.open_dataset(path + filename)
	pd = df.isel(time=slice(-28,None), lat=slice(-180,-45), lon=slice(-750,-500)) #last four weeks

	facet = pd[args.dataset].plot(x='lon', y='lat', col='time', col_wrap=7,robust=True,figsize=(11,8.5))
	facet.fig.set_size_inches( (16.5, 8.5) )	
	facet.fig.savefig(filename.replace('.nc',datetime.datetime.now().strftime('%b%d')+'.png'), dpi = (300))
	plt.close()