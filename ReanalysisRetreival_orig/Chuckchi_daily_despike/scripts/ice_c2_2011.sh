#!/bin/bash

# Purpose:
#       Script to run singlevarplot.py for each *.dat file in a directory

data_dir="/Users/bell/Data_Local/from_phyllis/iris/Chuckchi/ice_c2_2011/*.dat"
prog_dir="/Users/bell/Programs/Python/FOCI_Analysis/Chuckchi_daily_despike/"
out_dir="/Users/bell/Data_Local/from_phyllis/iris/Chuckchi/ice_c2_2011_images/"

for files in $data_dir
do
    names=(${files//\// })
	python ${prog_dir}singlevarplot.py ${files} ${out_dir} -date_fix 
done
