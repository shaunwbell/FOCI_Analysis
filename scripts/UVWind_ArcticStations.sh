#!/bin/bash

data_dir="/Users/bell/in_and_outbox/data_sets/reanalyis_data/NARR/daily/"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/"


#python NARR_daily_WindsSFCtemp_Station.py C1 70.836 163.124 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C2 71.230 164.105 2016 2019 --DataPath $data_dir
#python NARR_daily_WindsSFCtemp_Station.py C3 71.825 165.975 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C4 71.041 160.517 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C5 71.202 158.005 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C6 71.772 161.882 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C7 72.424 161.604 2010 2016
#python NARR_daily_WindsSFCtemp_Station.py C8 72.586 161.215 2010 2016
#python NARR_3hr_WindsSFCtemp_Station.py C9 72.467 156.550 2014 2017 --DataPath $data_dir
#python NARR_daily_WindsSFCtemp_Station.py C10 -
#python NARR_daily_WindsSFCtemp_Station.py C11 -
#python NARR_daily_WindsSFCtemp_Station.py C12 67.910 168.198 2010 2016
#python NARR_3hr_WindsSFCtemp_Station.py EARCTICSite 71.5 154.0 2014 2017 --DataPath $data_dir
#python NARR_3hr_WindsSFCtemp_Station.py WARCTICSite 75.0 170.0 2014 2017 --DataPath $data_dir


#anom dump
python NARR_daily_WindsSFCtemp_StationAnomdump.py C2 71.230 164.105 2016 2016 --DataPath $data_dir
