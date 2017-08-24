#!/usr/bin/env python

"""
 Background:
 --------
 oer_glider_dbload.py
 
 
 Purpose:
 --------
 load glider netcdf data into mysql database

 History:
 --------


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

data = pd.read_csv(args.sourcefile,delim_whitespace=True,skiprows=1,names=['obs','shipid'])

year,month = args.sourcefile.split('/')[-1].split('IDs_IMMA1_R3.0.0_')[-1].split('.txt')[0].split('-')

if pd.isnull(data['shipid']).any():
  print("{0} has shitty whitespace... skipping".format(args.sourcefile))
else:
  for index,row in data.iterrows():
      ICOADSdb.add_to_DB(table='ICOADS_ShipObs',Year=int(year),Month=int(month),ShipID=row.shipid,obs=row.obs)

ICOADSdb.close()