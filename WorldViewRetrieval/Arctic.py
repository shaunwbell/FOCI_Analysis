#!/usr/bin/env

"""
worldview_modis_wget.py

Purpose:
	Connect to ketch.pmel.noaa.gov (engineering run system) to retrieve prawler data from 2016 ITAE Bering Sea Mooring

"""
#System Stack
import datetime
import argparse

import wget

for i in range(0,91,1):
    dateutc=datetime.datetime(2019,5,10,0,0,0) + datetime.timedelta(days=i)
    datetimeutc_str=dateutc.strftime('%Y-%m-%dT%H:%M:%SZ')
    dateutc_str=dateutc.strftime('%Y-%m-%d')
    print(f'{dateutc_str}')
    url = f'https://wvs.earthdata.nasa.gov/api/v1/snapshot?REQUEST=GetSnapshot&TIME={datetimeutc_str}&BBOX=68.35436328070648,-169.4284675488015,74.44899877430296,-156.49032830625168&CRS=EPSG:4326&LAYERS=MODIS_Terra_CorrectedReflectance_TrueColor,Coastlines_15m,GHRSST_L4_MUR_Sea_Surface_Temperature,GHRSST_L4_MUR_Sea_Ice_Concentration&WRAP=day,x,day,day&FORMAT=image/tiff&WIDTH=2944&HEIGHT=1387&ts=1626400902553'

    wget.download(url, str(i).zfill(3)+'.tiff', bar=wget.bar_thermometer)