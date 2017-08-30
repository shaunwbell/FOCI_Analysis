# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 09:57:50 2016

@author: bell
"""

import datetime
import json
import csv
from natsort import natsorted

infile = '/Users/bell/Downloads/epickey.csv'
epic_file = csv.DictReader(open(infile,'rU'))

epic_keys = {}
for row in epic_file:
    epic_keys[row['KEY']] = row

print "{"
for ekey in natsorted(epic_keys.keys()):
    jd = json.dumps(epic_keys[ekey],sort_keys=True, indent=4)
    print "\"{0}\": {1},".format(epic_keys[ekey]['KEY'], jd)
print "}"    