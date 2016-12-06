#!/usr/bin/env

"""
    singlevarplot.py
    
    Plot single variable as a function of line number
    
"""
#System Stack
import datetime, sys
import argparse

#Science Stack
import numpy as np

# Visual Stack
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DateFormatter


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2015, 05, 13)
__modified__ = datetime.datetime(2015, 05, 13)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'plot', 'ice', 'one value'

"""------------------------General   Modules-------------------------------------------"""

parser = argparse.ArgumentParser(description='Plot single column ascii data')
parser.add_argument('DataPath', metavar='DataPath', type=str, help='full path to file')
parser.add_argument("OutputPath", metavar='OutputPath', type=str, help='full path to output file')
parser.add_argument("-date_fix",'--date_fix', action="store_true", help='assume file_ddMMMyyyy.dat and make sequential')
parser.add_argument("-plot",'--plot', action="store_true", help='save plots')

args = parser.parse_args()

print args.DataPath

################
### ingest file
data = {}
with open(args.DataPath) as f:

    for k, line in enumerate(f.readlines()):
        line = line.strip()

        data[k] = float(line)
        
#################
###Format date in file name
if args.date_fix:
    file_date = args.DataPath.split('/')[-1].split('_')[-1].split('.')[0]
    file_date = datetime.datetime.strptime(file_date,'%d%b%Y').strftime('%Y%m%d')
else:
    file_date = args.DataPath.split('/')[-1].split('_')[-1].split('.')[0]
    
################
#### Set output file
output_file = args.OutputPath + file_date + '_' +args.DataPath.split('/')[-1].split('.')[0] + '.png'

################
#### Print Derivative
deriv_temp = np.diff(data.values())
d_ind = np.where(np.abs(deriv_temp) > 5)
for dd in d_ind[0]:
    print dd+2, data[dd+1] #offset +1 for derivative, +1 for line number - value is just +1
    
################
#### Plot Derivative
deriv_temp = np.diff(data.values())
d_ind = np.where(np.abs(deriv_temp) < 5)
deriv_temp[d_ind] = 0

"""
ptitle = ("File: {0} \n ").format(args.DataPath.split('/')[-1])

fig = plt.figure(1)

ax2 = plt.subplot2grid((3, 1), (1, 0), colspan=1, rowspan=3)
p2 = ax2.plot(deriv_temp,'k-', markersize=2)
plt.ylabel('Values')
plt.xlabel('row #')

t = fig.suptitle(ptitle)
t.set_y(0.06)

fig.autofmt_xdate()
DefaultSize = fig.get_size_inches()
fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
plt.savefig(output_file, bbox_inches='tight', dpi = (100))

plt.close()
"""


    
################
#### Plot
if args.plot:
    ptitle = ("File: {0} \n ").format(args.DataPath.split('/')[-1])

    fig = plt.figure(1)

    ax2 = plt.subplot2grid((3, 1), (1, 0), colspan=1, rowspan=3)
    p2 = ax2.plot(data.keys(), data.values(),'k-', markersize=2)
    plt.ylabel('Values')
    plt.xlabel('row #')

    t = fig.suptitle(ptitle)
    t.set_y(0.06)

    fig.autofmt_xdate()
    DefaultSize = fig.get_size_inches()
    fig.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]) )
    plt.savefig(output_file, bbox_inches='tight', dpi = (100))

    plt.close()