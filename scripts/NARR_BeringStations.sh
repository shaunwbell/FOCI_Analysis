#!/bin/bash

data_dir="/Users/bell/in_and_outbox/data_sets/reanalyis_data/NARR/daily/"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/"


#python NARR_daily_WindsSFCtemp_Station.py M2 56.87 164.05 2005 2017 
python NARR_daily_WindsSFCtemp_Station.py M8 62.19 174.689 1980 2018 --DataPath ${data_dir}
