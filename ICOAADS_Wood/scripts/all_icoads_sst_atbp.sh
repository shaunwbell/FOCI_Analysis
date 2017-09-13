#Fall 2017 Deployment
data_dir="/Users/bell/in_and_outbox/2017/wood/icoads_bp_at_sst/*.csv"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/ICOAADS_Wood/"

for data_file in $data_dir
do
    echo ${data_file}
    python ${prog_dir}icoads_dbload_sst_atbp.py ${data_file}
done
