#!/bin/bash
# script to generate timeseries data from high resolution OI sst data for chosen sites

data_dir="/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/BOEM_arcwest_chaosx/NARR_SFCuv_daily/*C9*.nc"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/EcoFOCI_Utilities/"
outdir="/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/BOEM_arcwest_chaosx/NARR_SFCuv_daily/"

for files in $data_dir
do
    names=(${files//\// })
    outfile=${names[${#names[@]} - 1]}
    echo "processing file: $files"
	python ${prog_dir}nc2csv.py ${files} -timeseries -units_meta -EPIC WU_422 WV_423  >> ${outdir}${outfile}.csv
done