#Fall 2017 Deployment
data_dir="/Volumes/WDC_internal/Users/bell/in_and_outbox/2017/wood/ICOAADS_ship_obs/*.txt"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/ICOAADS_Wood/"

for data_file in $data_dir
do
    echo ${data_file}
    python ${prog_dir}icoads_dbload.py ${data_file}
done
