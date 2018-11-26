#!/bin/bash
# script to generate timeseries data from high resolution OI sst data for chosen sites

data_dir="/Users/bell/in_and_outbox/2018/stabeno/nov/m8_winds/*.nc"
prog_dir="/Users/bell/Programs/Python/EcoFOCI_Utilities/"
outdir="/Users/bell/in_and_outbox/2018/stabeno/nov/m8_winds/"

for files in $data_dir
do
    names=(${files//\// })
    outfile=${names[${#names[@]} - 1]}
    echo "processing file: $files"
	python ${prog_dir}nc2csv.py ${files} -timeseries -units_meta -EPIC WU_422 WV_423 AT_21 >> ${outdir}${outfile}.csv
done