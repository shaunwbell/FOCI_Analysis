# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 08:17:01 2016

@author: bell
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

### read in example dataset as a pandas dataframe.  it is timeseries data
# from the radiometer buoy (multiple timeseries values for shortwave radiation)
# as well as buoy orientation  
data = pd.read_csv('data/2016shiptest_OER_SWR_SPN1_tiltcorr.csv',skipinitialspace=True)
data['Date'] = pd.to_datetime(data['Date'],format='%Y-%m-%d %H:%M:%S')
#data['Date'] = pd.to_datetime(data['Date'],format='%m/%d/%Y %H:%M:%S')

styles = plt.style.available

"""
1 seaborn-darkgrid
seaborn-notebook
classic
seaborn-ticks
grayscale
bmh
seaborn-talk
dark_background
9 gplot
fivethirtyeight
seaborn-colorblind
seaborn-deep
seaborn-whitegrid
seaborn-bright
seaborn-poster
seaborn-muted
seaborn-paper
seaborn-white
seaborn-pastel
seaborn-dark
21 seaborn-dark-palette
"""
for s in styles:
    plt.style.use(s)
    data.hist()
    
## filter the data    
from scipy.fftpack import fft

N=69220
Y=fft(data[' corr_sza']-data[' sunzen'])
n=len(Y)
power = abs(Y[1:(n/2)])**2
nyquist=1./2
freq=np.array(range(n/2))/(n/2.0)*nyquist
period=1./freq