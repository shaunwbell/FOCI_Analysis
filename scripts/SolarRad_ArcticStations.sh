#!/bin/bash
# script to generate timeseries data from high resolution OI sst data for chosen sites

data_dir="/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/"

python NARR_rad_Station.py M2 56.864 164.0655 2017 2017
python NARR_rad_Station.py M8 62.2 174.678 2017 2017
python NARR_rad_Station.py CK 68.0 170.000 2017 2017

: '
python NARR_rad_Station.py C1 70.836 163.124 2010 2016
python NARR_rad_Station.py C2 71.230 164.105 2010 2016
python NARR_rad_Station.py C3 71.825 165.975 2010 2016
python NARR_rad_Station.py C4 71.041 160.517 2010 2016
python NARR_rad_Station.py C5 71.202 158.005 2010 2016
python NARR_rad_Station.py C6 71.772 161.882 2010 2016
python NARR_rad_Station.py C7 72.424 161.604 2010 2016
python NARR_rad_Station.py C8 72.586 161.215 2010 2016
python NARR_rad_Station.py C9 72.467 156.550 2010 2016
#python NARR_rad_Station.py C10 -
#python NARR_rad_Station.py C11 -
python NARR_rad_Station.py C12 67.910 168.198 2010 2016
'