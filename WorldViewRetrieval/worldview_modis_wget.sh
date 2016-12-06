#!/bin/bash

progdir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/WorldViewRetrieval/"

for i in `seq 174 194`;
do
        python ${progdir}EasternEqPacific.py ${i} VIIRS jpeg large
done  

img_dir="/Volumes/WDC_internal/Users/bell/Programs/Python/FOCI_Analysis/WorldViewRetrieval/*.jpeg"
for files in $img_dir
do
	convert ${files} -fill white -undercolor '#00000080' -pointsize 50 -gravity NorthEast -annotate +10+10 %t ${files}.jpg
done
