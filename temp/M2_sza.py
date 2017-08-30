#!/usr/bin/env python

"""
M2_solar_zenith_angle.py
 
"""

import numpy as np
import pandas as pandas
import general_utilities.solar_zenith_angle as sza

data = pandas.read_csv('/Users/bell/sites/bell/eFOCI_Mooring_logs/Mooring_CSV/15bs2a_ecf.csv',header=1,names=['time', 'fluor_3031', 'Fch_906'])
data['time'] = pandas.to_datetime(data['time'],format='%m/%d/%y %H:%M:%S')
sLength = len(data['time'])
data.loc[:,'sza'] = pandas.Series(np.random.randn(sLength), index=data.index)
data.loc[:,'az'] = pandas.Series(np.random.randn(sLength), index=data.index)

for index,row in data.iterrows():
    data['sza'][index], data['az'][index] = sza.solar_zenith(row['time'],lat='56.866667',lon='-164.05')
    if index % 100 == 0:
        print "{0} of {1}".format(index,sLength)
        
data.to_csv('/Users/bell/sites/bell/eFOCI_Mooring_logs/Mooring_CSV/15bs2a_ecf_sza.csv', columns=['time','fluor_3031','Fch_906','sza','az'], index=False)