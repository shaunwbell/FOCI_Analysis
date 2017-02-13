#!/bin/bash
# script to generate timeseries data from high resolution OI sst data for chosen sites

data_dir="/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/BOEM_arcwest_chaosx/NARR_solarrad_daily/*.nc"
prog_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/EcoFOCI_Utilities/"
outdir="/Volumes/WDC_internal/Users/bell/in_and_outbox/2016/BOEM_arcwest_chaosx/NARR_solarrad_daily/"

epic_vars="T_28 T2_35 S_41 S_42 ST_70 F_903 Trb_980 O_65 OST_62 CTDOXY_4221 CTDOST_4220 PAR_905 ATTN_55 Tr_904"

for files in $data_dir
do
    names=(${files//\// })
    outfile=${names[${#names[@]} - 1]}
    echo "processing file: $files"
	python ${prog_dir}nc2csv.py ${files} -timeseries -units_meta -EPIC Qs_133  >> ${outdir}${outfile}.csv
done