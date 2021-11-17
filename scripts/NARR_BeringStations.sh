#!/bin/bash

data_dir="/Volumes/MobileSSD/in_and_outbox/data_sets/reanalyis_data/NARR/3hr/"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/"


#v M2 56.87 164.05 2005 2017 
#python NARR_daily_WindsSFCtemp_Station.py M8 62.19 174.689 1980 2018 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py M8 62.19 174.689 1980 2019 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py M5 59.91 174.680 1980 2019 --DataPath ${data_dir}
#python NARR_3hr_WindsSFCtemp_Station.py M8 62.19 174.689 2018 2019 --DataPath ${data_dir} 

## Arctic
python NARR_daily_WindsSFCtemp_Station.py C2 71.230 164.105 2015 2020 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py C2 71.230 164.105 2018 2019 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py C9 72.467 156.550  2019 2019 --DataPath ${data_dir}

### GOA Stations
#python NARR_3hr_WindsSFCtemp_Station.py Shumagin 55.0 160.5 2015 2017 --DataPath ${data_dir} -sm


### 
#python NARR_daily_WindsSFCtemp_Station.py IF 58.669 168.319 1997 2000 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py M2 56.867 164.05 1997 2000 --DataPath ${data_dir}
#python NARR_daily_WindsSFCtemp_Station.py KD 59.063 164.173 1997 2000 --DataPath ${data_dir}
