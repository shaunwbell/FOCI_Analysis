Readme
======

Purpose
-------

The routines in this directory have one of the following purposes:   

* average mooring data to 24hr bins and compare to NARR/NCEP V2 24hr data (nearest point)   
* average mooring data to 12hr bins and average 3hr/6hr NARR/NCEP data to 12hr bins to compare   
* average mooring data to 3hr bins and compare to NARR 3hr data (nearest point)   
* prepare EPIC flavored NetCDF files with u, v winds from NARR 3hr data (smoothed with 3point triangular filter)   

### Outputs   

The outputs fall into one of the following categories:   

* 1:1 comparison plots of moorings with reanalysis data for u,v components and along/across shore components   
* Plot topography/bathymetry of region of interest   
* Plot comparisons of NARR and NCEP analysis for same time and place (wind field)   
* Timeseries of u,v components for all three data sources   
* an attempt at speed, direction analysis as well (less useful and not being updated)   

### Other Notes:   

* Mooring data that spans multiple seasons can be read in as one file (MFDataset)
* Nearest point uses haversine formula   
* EPIC NetCDF generation is a useful class routine definition and should be expanded on

### References:   

* EcoFOCI - Mooring data   
* NCEP - NARR   
* NCEP - Reanalysis V2
I