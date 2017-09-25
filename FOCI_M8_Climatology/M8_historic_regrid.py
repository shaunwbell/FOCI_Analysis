#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:45:41 2017

@author: bell
"""

import numpy as np
import pandas as pd

data = pd.read_excel('/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/stabeno/BeringM8_Historic/M8_50kmrad_allhistoric_profiles_populated.xlsx',1)

dgrp = data.groupby('CastID')

for d in dgrp.groups.keys():
    td =dgrp.get_group(d)

    nd = np.arange(0,np.ceil(td['Depth [m]'].max()))
    nt = np.interp(np.arange(0,np.ceil(td['Depth [m]'].max())),td['Depth [m]'].values,td['Temperature [degrees_C]'].values)
    ns = np.interp(np.arange(0,np.ceil(td['Depth [m]'].max())),td['Depth [m]'].values,td['Salinity [psu]'].values)
    
    for i,v in enumerate(nd):
        print d, td['yyyy-mm-ddThh:mm:ss.sss'][td['yyyy-mm-ddThh:mm:ss.sss'].keys()[0]],td['Latitude [degrees_north]'][td['Latitude [degrees_north]'].keys()[0]],td['Longitude [degrees_east]'][td['Longitude [degrees_east]'].keys()[0]],v, nt[i], ns[i]