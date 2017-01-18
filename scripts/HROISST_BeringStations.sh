#!/bin/bash
# script to generate timeseries data from high resolution OI sst data for chosen sites

data_dir="/Volumes/WDC_internal/Users/bell/Data_Local/sst/NOAA_OI_SST_V2/"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/"


python HROISST_Station.py M4 57.8611 168.884 1981 2016 -scf -sep

python HROISST_Station.py M8 62.19 174.689 1981 2016 -scf -sep

python HROISST_Station.py M2 56.87 164.05 1981 1995 -scf -sep