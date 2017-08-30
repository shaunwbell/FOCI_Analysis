#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 10:37:24 2017

@author: bell
"""

import pandas as pd

f='/Volumes/WDC_internal/Users/bell/Downloads/AW_and_CHX_GAM_AllComponents.xlsx'

data = pd.read_excel(f,sheetname='C9', header=3)

column_to_interp=['F','T','S','U','V','W','X','R','u','t','n']

for var in column_to_interp:
    data[var].interpolate(limit=7,inplace=True)

print data.to_csv()