# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:42:45 2017

@author: bell
"""

import xarray as xr
import glob
import argparse

"""----------------------------- Main -------------------------------------"""

# parse incoming command line options
parser = argparse.ArgumentParser(description='Map')
parser.add_argument('sourcedir', metavar='sourcedir', type=str, help='full path to file')
args = parser.parse_args()

path = args.sourcedir
files = glob.glob(path + '*cf.nc')

for yearf in files:
	ds = xr.open_dataset(yearf)
	dsr = ds.resample('3MS',dim='time', how='mean')
	print dsr.to_dataframe().to_csv()
	ds.close()
	

"""
ds = xr.open_mfdataset(files)
dsr = ds.resample('MS',dim='time', how='mean')
print dsr.T_20.to_dataframe()
ds.close()
"""