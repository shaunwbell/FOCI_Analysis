#!/bin/bash

data_dir="/Volumes/WDC_internal/Users/bell/Data_Local/Reanalysis_Files/NARR/daily/"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/"


#v M2 56.87 164.05 2005 2017 
#python NARR_daily_WindsSFCtemp_Station.py M8 62.19 174.689 1980 2018 --DataPath ${data_dir}
#python NARR_daily_WindsTemp_atheight.py M8 62.19 174.689 1980 2018 --DataPath ${data_dir}
python NARR_daily_WindsSFCtemp_Station.py C2 71.230 164.105 2014 2017 --DataPath ${data_dir}