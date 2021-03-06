#!/usr/bin/env

"""
 Background:
 --------
 PlotMooringMultiInst.py
 
 
 Purpose:
 --------
 Plot Multiple Timeseries on same panel.  MULTIPLOT Overlay

 	Use Case1:
 		Plot the same parameter (eg temperature) from multiple instruments on same mooring from 
 		different depths. e.g. temperature from every MTR depth of a mooring deployment

 	Use Case2:
 		Plot the same parameter from different moorings at the same depth (not restricted to this).
 		e.g. temperature from the instrument closest to the surface over multiple deployments (onrunning SST plots)

 	Use Case3 (with -ctd flag):
 		Plot the discrete point from a ctd cast (nearby is most relevant) for QC puposes

 Modifications:
 --------------

 2016-09-16: SW Bell - Add support for parsing yaml files and translating between yaml and json/pyini
 					Begin code cleanup from previous iterations of the routine.  Merge so that one program can provide ctd cal
 					overlays.
  

"""

#System Stack
import datetime, sys, os
import argparse

#Science Stack
from netCDF4 import Dataset
import numpy as np

# Visual Stack
import matplotlib as mpl
mpl.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, WeekdayLocator, MonthLocator, DayLocator, HourLocator, DateFormatter
import matplotlib.ticker as ticker

# User Stack
from io_utils import ConfigParserLocal
from calc.EPIC2Datetime import EPIC2Datetime, get_UDUNITS
from io_utils.EcoFOCI_netCDF_read import EcoFOCI_netCDF

__author__   = 'Shaun Bell'
__email__	= 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2014, 9, 11)
__modified__ = datetime.datetime(2016, 9, 5)
__version__  = "0.1.0"
__status__   = "Development"
__keywords__ = 'Mooring', 'comparisons', 'Cruise', 'plots'


"""--------------------------------main Routines---------------------------------------"""

parser = argparse.ArgumentParser(description='SBE56 plotting')
parser.add_argument('PointerFile', 
	metavar='PointerFile', 
	type=str, 
	help='full path to pointer file')
parser.add_argument("-mt",'--manual_timebounds', 
	nargs='+', 
	type=str, 
	help='set times to specified values (d-m-Y)')
parser.add_argument("-md",'--manual_databounds', 
	nargs='+', 
	type=float, 
	help='set databounds to specified values')
parser.add_argument("-multi",'--multiplot_overlay', 
	action="store_true", 
	help='plot multiple mooring data on one panel')
parser.add_argument("-ctd",'--ctd_calibration_plots', 
	action="store_true", 
	help='plot CTD calibration point on timeseries')

args = parser.parse_args()

"""---------------------------------------------------------------------------------------
Get parameters from specified pointerfile - 
an example is shown in the header description of
this program.  It can be of the .pyini (json) form or .yaml form

"""
if args.PointerFile.split('.')[-1] == 'pyini':
	pointer_file = ConfigParserLocal.get_config(args.PointerFile)
elif args.PointerFile.split('.')[-1] == 'yaml':
	pointer_file = ConfigParserLocal.get_config_yaml(args.PointerFile)
else:
	print "PointerFile format not recognized"
	sys.exit()

MooringDataPath = pointer_file['mooring_data_path']
MooringID = pointer_file['MooringID']
color_options = pointer_file['colors']
label = pointer_file['legend']
legend_loc = pointer_file['legend_loc']
legend_off = pointer_file['legend_off']
datatype = pointer_file['dtype']
plot_var = pointer_file['EPIC_Key']
plot_var_ctd = pointer_file['EPIC_Key_ctd']
LocatorInterval = pointer_file['Date_Ticks']
Ylabel = pointer_file['Ylabel']
output_type = pointer_file['output_type']

MooringDataPath = pointer_file['mooring_data_path']
files = pointer_file['mooring_files']
files_path = [a+b for a,b in zip(MooringDataPath,files)]

CTDDataPath = pointer_file['ctd_data_path']
ctd_files = pointer_file['ctd_files']
ctd_files_path = [a+b for a,b in zip(CTDDataPath,ctd_files)]

### some mpl specif settings for fonts and plot style
mpl.rcParams['svg.fonttype'] = 'none'
plt.style.use(pointer_file['plot_stylesheet'])
#seaborn-poster -- fonts are smaller
#ggplot -- grey border, better axis frame
#bmh -- slightly heavier than ggplot for line weights


"""---------------------------------------------------------------------------------------
		Plot Multiple Mooring Datastreams on one panel
"""
databounds={}

if args.multiplot_overlay:
	### set title for plot
	ptitle = ("Plotted on: {timestr} \n from {mooringid} ").format(timestr=datetime.datetime.now().strftime('%Y/%m/%d %H:%M'), 
													 mooringid=MooringID )

	### initialize plot
	fig = plt.figure()
	plt.subplot2grid((3, 1), (1, 0), colspan=1, rowspan=3)

	### set arbitrary max and min bounds to be changed later based on data bounds
	databounds['max_t'] = 0
	databounds['min_t'] = 100000000
	databounds['max_v'] = -50
	databounds['min_v'] = 50
	label_thin = []

	### cycle through all files, retrieve data and plot
	print files_path
	for ind, ncfile in enumerate(files_path):
		print "Working on {activefile}".format(activefile=ncfile)

		#open/read netcdf files
		df = EcoFOCI_netCDF(ncfile)
		global_atts = df.get_global_atts()
		vars_dic = df.get_vars()
		ncdata = df.ncreadfile_dic()
		df.close()

		nctime = get_UDUNITS(EPIC2Datetime(ncdata['time'],ncdata['time2']),'days since 0001-01-01') + 1.

		#find and replace missing values with nans so they don't plot
		for var in plot_var:
			try:
				ncdata[var][np.where(ncdata[var] >1e30)] = np.nan
				label_thin = label_thin + [label[ind]]
			except KeyError:
				pass

		#Plot data
		plt.hold(True)
		for var in plot_var:
			if var in ['T_25']:
				try:
					ncdata['ICEC_2025'][np.where(ncdata['ICEC_2025'] >1e30)] = np.nan
					ncdata[var][np.where(ncdata['ICEC_2025'] > 1)] = np.nan
					plt.plot(nctime, ncdata[var][:,0,0,0],color_options[ind],linewidth=0.25)
				except KeyError: #if the file doesn't have the specified epic_key it will through an exception
					print "Failed to plot {0}".format(var)
					continue
			else:
				try:
					plt.plot(nctime, ncdata[var][:,0,0,0],color_options[ind],linewidth=0.25)
				except KeyError: #if the file doesn't have the specified epic_key it will through an exception
					print "Failed to plot {0}".format(var)
					continue
		#setup bouds
		for var in plot_var:
			try:
				if nctime.max() > databounds['max_t']:
					databounds['max_t'] = nctime.max()
				if nctime.min() < databounds['min_t']:
					databounds['min_t'] = nctime.min()
				if np.nanmax(ncdata[var][:,0,0,0]) > databounds['max_v']:
					databounds['max_v'] = np.nanmax(ncdata[var][:,0,0,0])
				if np.nanmin(ncdata[var][:,0,0,0]) < databounds['min_v']:
					databounds['min_v'] = np.nanmin(ncdata[var][:,0,0,0])
			except KeyError:
				pass

	#set bounds if estabilshed by user
	if args.manual_timebounds:
		databounds['min_t'] = datetime.datetime.strptime(args.manual_timebounds[0],'%Y-%m-%d').toordinal()
		databounds['max_t'] = datetime.datetime.strptime(args.manual_timebounds[1],'%Y-%m-%d').toordinal()

	#set bounds if estabilshed by user
	if args.manual_databounds:
		databounds['min_v'] = args.manual_databounds[0]
		databounds['max_v'] = args.manual_databounds[1]
		

	ax2 = plt.gca()
	ax2.set_ylim(databounds['min_v'],databounds['max_v'])
	ax2.set_xlim([databounds['min_t'],databounds['max_t']])

	if not legend_off:
		leg = ax2.legend(label_thin, loc=legend_loc, ncol=6, fontsize=8)
		for legobj in leg.legendHandles:
			legobj.set_linewidth(2.0)
		
	plt.ylabel(Ylabel)
	if LocatorInterval == 'multi_year':
		ax2.xaxis.set_major_locator(YearLocator())
		ax2.xaxis.set_minor_locator(MonthLocator(bymonth=6))
		ax2.xaxis.set_major_formatter(ticker.NullFormatter())
		ax2.xaxis.set_minor_formatter(DateFormatter('%Y'))
		ax2.tick_params(axis='both', which='minor', labelsize=12)
	else:
		ax2.xaxis.set_major_locator(MonthLocator())
		ax2.xaxis.set_minor_locator(MonthLocator(bymonth=[1,3,5,7,9,11], bymonthday=15))
		ax2.xaxis.set_major_formatter(ticker.NullFormatter())
		ax2.xaxis.set_minor_formatter(DateFormatter('%b %y'))
		ax2.tick_params(axis='both', which='minor', labelsize=12)

	t = fig.suptitle(ptitle, fontsize=8)
	t.set_y(0.03)

	#fig.autofmt_xdate()
	DefaultSize = fig.get_size_inches()
	fig.set_size_inches( (DefaultSize[0], DefaultSize[1]) )
	plt.savefig('images/'+ MooringID + '_'+plot_var[0]+'_'+datatype+'.'+output_type, bbox_inches='tight', dpi = (300))
	plt.close()

