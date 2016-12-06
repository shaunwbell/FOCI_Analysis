# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 06:42:17 2016

@author: bell
"""

import datetime
import numpy as np

import gsw

from io_utils.EcoFOCI_netCDF_read import EcoFOCI_netCDF

#godas is a function of time, depth, lat, lon (12, 40, 418, 360)

godas_ptemp = '/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/GODAS/pottmp.1980.nc'
godas_sal = '/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/GODAS/salt.1980.nc'


#GODAS PTEMP
df = EcoFOCI_netCDF(godas_ptemp)
global_atts = df.get_global_atts()
vars_dic = df.get_vars()
gd_ptmp = df.ncreadfile_dic()
df.close()

#GODAS SAL
df = EcoFOCI_netCDF(godas_sal)
global_atts = df.get_global_atts()
vars_dic = df.get_vars()
gd_sal = df.ncreadfile_dic()
df.close()


#ABS Sal f(sal, pres, lat, lon)
#pressure needs to be determined from depth

ABS_SA = np.ones_like(gd_sal['salt'])
CONS_T = np.ones_like(gd_sal['salt'])
for time_ind,month in enumerate(gd_sal['time']):
	print time_ind
	for level_ind, level in enumerate(gd_sal['level']):
		print level
		for lat_ind,lat in enumerate(gd_sal['lat']):
			if lat_ind % 100 == 0:
				print "{0} of {1}".format(lat_ind, ABS_SA.shape[2])
			press = gsw.p_from_z(-1.0*level,lat)
			for lon_ind,lon in enumerate(gd_sal['lon']):
				if lon_ind % 100 == 0:
					print "{0} of {1}".format(lon_ind, ABS_SA.shape[3])
				ABS_SA[time_ind, level_ind, lat_ind, lon_ind] = \
					gsw.SA_from_SP(gd_sal['salt'][time_ind,level_ind,
									lat_ind,lon_ind], press,lat,lon)
				CONS_T[time_ind, level_ind, lat_ind, lon_ind] = \
					gsw.CT_from_pt(ABS_SA[time_ind, level_ind, lat_ind, lon_ind], 
									gd_ptmp['pottmp'][time_ind,level_ind,
									lat_ind,lon_ind])


