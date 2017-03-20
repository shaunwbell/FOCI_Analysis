#!/usr/bin/env python
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
import cmocean

__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2017, 1, 30)
__modified__ = datetime.datetime(2017, 1, 30)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'xarray','SST','anom','NOAA HRES OI V2 SST'

"""--------------------------------helper Routines---------------------------------------"""


"""--------------------------------main Routines---------------------------------------"""

parser = argparse.ArgumentParser(description='NASA WorldView MODIS Retrieval')
parser.add_argument('year', metavar='year', type=str, help='year of data set')               
parser.add_argument('dataset', metavar='dataset', type=str, help='dataset ("anom","err","sst, icec")' )
parser.add_argument('-plot','--plot', type=str, help='path to save plot of last 4 weeks to' )
                  
args = parser.parse_args()

if args.dataset == 'sst':
	var = 'sst'
	cmap = cmocean.cm.thermal
	param = 'mean'
	plotvar = 'sst'
elif args.dataset == 'anom':
	var = 'sst'
	cmap = 'RdBu_r'
	param = args.dataset
	plotvar = 'anom'
elif args.dataset == 'icec':
	var = 'icec'
	cmap = 'RdBu'
	param = 'mean'
	plotvar = 'icec'
else:
	var = args.dataset
	cmap = 'RdBu_r'
	param = args.dataset
	plotvar = args.dataset

ftp = FTP('ftp.cdc.noaa.gov',user='anonymous')
ftp.cwd('Datasets/noaa.oisst.v2.highres')
path = '/Volumes/WDC_internal/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/'
filename = var+'.day.'+param+'.'+args.year+'.v2.nc'
# download
fhandle = open(path + filename, 'wb')
ftp.retrbinary('RETR ' + filename, fhandle.write)                                       # Imprimimos por pantalla lo que estamos descargando        #fhandle.close()
fhandle.close()  
ftp.close()

if args.plot:
	df = xa.open_dataset(path + filename)
	pd = df.isel(time=slice(-28,None), lat=slice(-180,-45), lon=slice(-750,-500)) #last four weeks

	facet = pd[plotvar].plot(x='lon', y='lat', col='time', col_wrap=7,robust=True,figsize=(11,8.5), cmap=cmap)
	facet.fig.set_size_inches( (16.5, 8.5) )	
	facet.fig.savefig(args.plot + filename.replace('.nc',datetime.datetime.now().strftime('%b%d')+'.png'), dpi = (300))
	plt.close()