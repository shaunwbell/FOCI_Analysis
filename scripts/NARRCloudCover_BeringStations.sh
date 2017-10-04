#!/bin/bash

data_dir="/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/3hourly/"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/"


#python NARR_3hr_WindsSFCtemp_Station.py M8 62.19 174.689 2016 2016 -sm

python NARR_daily_3CloudCover_Station.py M2 56.87 164.05 2005 2016 
python NARR_daily_3CloudCover_Station.py M8 62.19 174.689 2005 2016 
