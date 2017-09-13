#!/usr/bin/env python

"""
 Background:
 --------
 icoads_dbload_sst.py
 
 
 Purpose:
 --------
 load icoaads data into mysql database

 History:
 --------


--
YR,MO,DY,HR,LAT,LON,ID,SST,DCK,SID,PT,UID,RN1,RN2,RN3,RSA,IRF


"""
import argparse, os

from io_utils import ConfigParserLocal
from io_utils.EcoFOCI_db_io import EcoFOCI_db
import pandas as pd

"""-------------------------------- Main -----------------------------------------------"""

parser = argparse.ArgumentParser(description='Load ICOADS Text File')
parser.add_argument('sourcefile', metavar='sourcefile', type=str,
               help='complete path to icoads file')
args = parser.parse_args()



###
#
# load database
config_file = 'EcoFOCI_config/db_config/db_config_sideprojects.pyini'
ICOADSdb = EcoFOCI_db()
(db,cursor) = ICOADSdb.connect_to_DB(db_config_file=config_file)

data = pd.read_csv(args.sourcefile, error_bad_lines=False)
data = data.where((pd.notnull(data)), None)

for index,row in data.iterrows():
  ICOADSdb.add_to_DB(table='ICOADS_SST_ATBP',YR=row.YR,MO=row.MO,DY=row.DY,HR=row.HR,LAT=row.LAT,LON=row.LON,
  						ID=row.ID,AT=row.AT,SLP=row.SLP)
  if index % 1000 == 0:
  	print("Inserted row {index}".format(index=index))

ICOADSdb.close()