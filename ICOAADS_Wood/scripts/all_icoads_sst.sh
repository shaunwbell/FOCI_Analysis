#Fall 2017 Deployment
data_dir="/Users/bell/in_and_outbox/2017/wood/icoaads_all/*.csv"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/ICOAADS_Wood/"

for data_file in $data_dir
do
    echo ${data_file}
    python ${prog_dir}icoads_dbload_sst.py ${data_file}
done
