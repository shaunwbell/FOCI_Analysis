#!/usr/bin/env 
"""
 ncutilities.py
 

   Using Anaconda packaged Python
"""
import datetime
from netCDF4 import Dataset 
from netCDF4 import MFDataset 


__author__   = 'Shaun Bell'
__email__    = 'shaun.bell@noaa.gov'
__created__  = datetime.datetime(2013, 12, 20)
__modified__ = datetime.datetime(2014, 03, 14) #<-- PI DAY :)
__version__  = "0.1.1"
__status__   = "Development"

"""------------------------------NC open/read/close------------------------------------"""

def ncopen(ncfile):
    """
    Parameters
    ----------
    TODO
    
    Returns
    -------
    TODO
              
    """
    nchandle = Dataset(ncfile,'r')
    return nchandle
        
def ncclose(nchandle):
    """
    Parameters
    ----------
    TODO
    
    Returns
    -------
    TODO
              
    """
    nchandle.close()

def get_global_atts(nchandle):
    """returns global attribute names"""
    att_names = nchandle.ncattrs()
    return att_names

def get_vars(nchandle):
    """returns just the variables"""
    return nchandle.variables.keys()
    
def get_geocoords(nchandle, lat='lat', lon='lon'):
    
    data = {}
    
    for j, v in enumerate([lat, lon]):
        data[v] = nchandle.variables[v][:]
        
    return (data)
    
def ncreadfile_dic(nchandle, params):
    """returns all data for all variables"""
    data = {}
    for j, v in enumerate(params): 
        if v in nchandle.variables.keys(): #check for nc variable
                data[v] = nchandle.variables[v][:]
        else: #if parameter doesn't exist fill the array with zeros
            data[v] = None
    return (data)
    
def ncreadfile_dic_slice(nchandle, params, height_ind=None, lat_ind=None, lon_ind=None):
    """returns slice of data for all times but for specified height/lat/lon indicies"""
    data = {}
    if height_ind == None:
        for j, v in enumerate(params): 
            try: #check for nc variable
                    data[v] = nchandle.variables[v][:,lat_ind,lon_ind]

            except ValueError: #if parameter is not of expected dimensions
                data[v] = nchandle.variables[v][:]
    else:
        for j, v in enumerate(params): 
            try: #check for nc variable
                    data[v] = nchandle.variables[v][:,:,lat_ind,lon_ind]

            except ValueError: #if parameter is not of expected dimensions
                data[v] = nchandle.variables[v][:]

    return data    
    
    
"""----------------------------------NC Save-------------------------------------------"""

"""-------------------------------EPIC Standard----------------------------------------"""

class EPIC_NC(object):
    """ Class instance to generate a NetCDF file.  
    

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = EPIC_NC()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'
    def __init__(self, savefile='data/test.nc', data=None):
        
        self.savefile = savefile
    
    def file_create(self):
        """ sets up filename, read/write permissions, netcdf format """
        
        rootgrpID = Dataset(self.savefile, EPIC_NC.nc_read, format=EPIC_NC.nc_format)
        self.rootgrpID = rootgrpID
        return ( rootgrpID )
        
    def sbeglobal_atts(self, DATA_CMNT='', coord_system="GEOGRAPHICAL", Water_Mass="G"):
        """
        Assumptions
        -----------
        
        Follows Mooring Format
        
        seabird related global attributes found in DataFrame.header list
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.DATA_CMNT = DATA_CMNT
        self.rootgrpID.COORD_SYSTEM = coord_system
        self.rootgrpID.WATER_MASS = Water_Mass
        
    def PMELglobal_atts(self, Barometer=9999, Wind_Dir=999, Wind_Speed=99,
                        Air_Temp=99.9, Water_Depth=9999, file_name='', Prog_Cmnt='', Edit_Cmnt='', Station_Name=''):
        """
        Assumptions
        -----------
        
        Format of DataFrame.name = 'dy1309l1_ctd001'
        
        seabird related global attributes found in DataFrame.header list
        
        Options
        -------
        
        Todo
        -----
        
        Retrieve PMEL header information from '@' comments in .cnv file or from separate
        header.txt file
        """        
        
        #From PMELheader
        
        self.rootgrpID.BAROMETER = Barometer
        self.rootgrpID.WIND_DIR = Wind_Dir
        self.rootgrpID.WIND_SPEED = Wind_Speed
        self.rootgrpID.AIR_TEMP = Air_Temp
        self.rootgrpID.WATER_DEPTH = Water_Depth
        self.rootgrpID.STATION_NAME = Station_Name
        self.rootgrpID.EPIC_FILE_GENERATOR = file_name
        self.rootgrpID.PROG_CMNT = Prog_Cmnt
        self.rootgrpID.EDIT_CMNT = Edit_Cmnt
        
        pass
        
    def dimension_init(self, len_time=1, len_depth=1, len_lat=1, len_lon=1 ):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time', 'depth', 'latitude', 'longitude']
        
        self.rootgrpID.createDimension( self.dim_vars[0], len_time ) #time
        self.rootgrpID.createDimension( self.dim_vars[1], len_depth ) #depth
        self.rootgrpID.createDimension( self.dim_vars[2], len_lat ) #lat
        self.rootgrpID.createDimension( self.dim_vars[3], len_lon ) #lon
        
        
    def variable_init(self):
        """
        U, V components of wind (m/s) - 'WU_422', 'WV_423'
        """

        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = ['WU_422', 'WV_423','AT_21'], ['WU', 'WV','AT'], ["WIND U (M/S)             ","WIND V (M/S)             ","AIR TEMPERATURE (C)"]
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = ['u','v','atemp'], ['','',''], ['m s-1','m s-1','C'], [422,423,21]
         


        rec_vars = ['time','time2','depth','latitude','longitude'] + rec_vars

        rec_var_name = ['', '', '', '', ''] + rec_var_name
        rec_var_longname = ['', '', '', '', ''] + rec_var_longname
        rec_var_generic_name = ['', '', '', '', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['f10.0', 'f10.0', 'f10.1', 'f10.4', 'f10.4'] + rec_var_FORTRAN
        rec_var_units = ['True Julian Day', 'msec since 0:00 GMT','m','degree_north','degree_west'] + rec_var_units
        rec_var_type= ['i4', 'i4'] + ['f4' for spot in rec_vars[2:]]
        rec_var_strtype= ['EVEN', 'EVEN', 'EVEN', 'EVEN', 'EVEN'] + ['' for spot in rec_vars[5:]]
        rec_epic_code = [624, 624,1,500,501] + rec_var_epic
        
        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[0]))#time2
        var_class.append(self.rootgrpID.createVariable(rec_vars[2], rec_var_type[2], self.dim_vars[1]))#depth
        var_class.append(self.rootgrpID.createVariable(rec_vars[3], rec_var_type[3], self.dim_vars[2]))#lat
        var_class.append(self.rootgrpID.createVariable(rec_vars[4], rec_var_type[4], self.dim_vars[3]))#lon
        
        for i, v in enumerate(rec_vars[5:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+5], rec_var_type[i+5], self.dim_vars))

        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, depth_level=-10., latitude=None, longitude=None, time1=None, time2=None, CastLog=False):
        """ """
        self.var_class[0][:] = time1
        self.var_class[1][:] = time2

        self.var_class[2][:] = depth_level #self.data[pressure_var].values
        self.var_class[3][:] = latitude
        self.var_class[4][:] = longitude #PMEL standard direction W is +

    def add_data(self, parameter, value):
        """ """

        di = self.rec_vars.index(parameter)
        self.var_class[di][:] = value
            
        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.History = self.History + ' ' + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                    + ' ' + new_history + '\n'
                    
    def close(self):
        self.rootgrpID.close()
    
class EPIC_NC_rad(object):
    """ Class instance to generate a NetCDF file.  
    

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = EPIC_NC_rad()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'
    def __init__(self, savefile='data/test.nc', data=None):
        
        self.savefile = savefile
    
    def file_create(self):
        """ sets up filename, read/write permissions, netcdf format """
        
        rootgrpID = Dataset(self.savefile, EPIC_NC_rad.nc_read, format=EPIC_NC_rad.nc_format)
        self.rootgrpID = rootgrpID
        return ( rootgrpID )
        
    def sbeglobal_atts(self, DATA_CMNT='', coord_system="GEOGRAPHICAL", Water_Mass="G"):
        """
        Assumptions
        -----------
        
        Follows Mooring Format
        
        seabird related global attributes found in DataFrame.header list
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.DATA_CMNT = DATA_CMNT
        self.rootgrpID.COORD_SYSTEM = coord_system
        self.rootgrpID.WATER_MASS = Water_Mass
        
    def PMELglobal_atts(self, Barometer=9999, Wind_Dir=999, Wind_Speed=99,
                        Air_Temp=99.9, Water_Depth=9999, file_name='', Prog_Cmnt='', Edit_Cmnt='', Station_Name=''):
        """
        Assumptions
        -----------
        
        Format of DataFrame.name = 'dy1309l1_ctd001'
        
        seabird related global attributes found in DataFrame.header list
        
        Options
        -------
        
        Todo
        -----
        
        Retrieve PMEL header information from '@' comments in .cnv file or from separate
        header.txt file
        """        
        
        #From PMELheader
        
        self.rootgrpID.BAROMETER = Barometer
        self.rootgrpID.WIND_DIR = Wind_Dir
        self.rootgrpID.WIND_SPEED = Wind_Speed
        self.rootgrpID.AIR_TEMP = Air_Temp
        self.rootgrpID.WATER_DEPTH = Water_Depth
        self.rootgrpID.STATION_NAME = Station_Name
        self.rootgrpID.EPIC_FILE_GENERATOR = file_name
        self.rootgrpID.PROG_CMNT = Prog_Cmnt
        self.rootgrpID.EDIT_CMNT = Edit_Cmnt
        
        pass
        
    def dimension_init(self, len_time=1, len_depth=1, len_lat=1, len_lon=1 ):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time', 'depth', 'latitude', 'longitude']
        
        self.rootgrpID.createDimension( self.dim_vars[0], len_time ) #time
        self.rootgrpID.createDimension( self.dim_vars[1], len_depth ) #depth
        self.rootgrpID.createDimension( self.dim_vars[2], len_lat ) #lat
        self.rootgrpID.createDimension( self.dim_vars[3], len_lon ) #lon
        
        
    def variable_init(self):
        """
        U, V components of wind (m/s) - 'WU_422', 'WV_423'
        """

        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = ['Qs_133'], ['QS'], ["SHORTWAVE RADIATION"]
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = [''], [''], ['w m-2'], [133]
         


        rec_vars = ['time','time2','depth','latitude','longitude'] + rec_vars

        rec_var_name = ['', '', '', '', ''] + rec_var_name
        rec_var_longname = ['', '', '', '', ''] + rec_var_longname
        rec_var_generic_name = ['', '', '', '', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['f10.0', 'f10.0', 'f10.1', 'f10.4', 'f10.4'] + rec_var_FORTRAN
        rec_var_units = ['True Julian Day', 'msec since 0:00 GMT','m','degree_north','degree_west'] + rec_var_units
        rec_var_type= ['i4', 'i4'] + ['f4' for spot in rec_vars[2:]]
        rec_var_strtype= ['EVEN', 'EVEN', 'EVEN', 'EVEN', 'EVEN'] + ['' for spot in rec_vars[5:]]
        rec_epic_code = [624, 624,1,500,501] + rec_var_epic
        
        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[0]))#time2
        var_class.append(self.rootgrpID.createVariable(rec_vars[2], rec_var_type[2], self.dim_vars[1]))#depth
        var_class.append(self.rootgrpID.createVariable(rec_vars[3], rec_var_type[3], self.dim_vars[2]))#lat
        var_class.append(self.rootgrpID.createVariable(rec_vars[4], rec_var_type[4], self.dim_vars[3]))#lon
        
        for i, v in enumerate(rec_vars[5:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+5], rec_var_type[i+5], self.dim_vars))

        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, depth_level=-10., latitude=None, longitude=None, time1=None, time2=None, CastLog=False):
        """ """
        self.var_class[0][:] = time1
        self.var_class[1][:] = time2

        self.var_class[2][:] = depth_level #self.data[pressure_var].values
        self.var_class[3][:] = latitude
        self.var_class[4][:] = longitude #PMEL standard direction W is +

    def add_data(self, parameter, value):
        """ """

        di = self.rec_vars.index(parameter)
        self.var_class[di][:] = value
            
        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.History = self.History + ' ' + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                    + ' ' + new_history + '\n'
                    
    def close(self):
        self.rootgrpID.close()
    
class EPIC_NC_SST(object):
    """ Class instance to generate a NetCDF file.  
    

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = EPIC_NC()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'
    def __init__(self, savefile='data/test.nc', data=None):
        
        self.savefile = savefile
    
    def file_create(self):
        """ sets up filename, read/write permissions, netcdf format """
        
        rootgrpID = Dataset(self.savefile, EPIC_NC.nc_read, format=EPIC_NC.nc_format)
        self.rootgrpID = rootgrpID
        return ( rootgrpID )
        
    def sbeglobal_atts(self, DATA_CMNT='', coord_system="GEOGRAPHICAL", Water_Mass="G"):
        """
        Assumptions
        -----------
        
        Follows Mooring Format
        
        seabird related global attributes found in DataFrame.header list
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.DATA_CMNT = DATA_CMNT
        self.rootgrpID.COORD_SYSTEM = coord_system
        self.rootgrpID.WATER_MASS = Water_Mass
        
    def PMELglobal_atts(self, Barometer=9999, Wind_Dir=999, Wind_Speed=99,
                        Air_Temp=99.9, Water_Depth=9999, file_name='', Prog_Cmnt='', Edit_Cmnt='', Station_Name=''):
        """
        Assumptions
        -----------
        
        Format of DataFrame.name = 'dy1309l1_ctd001'
        
        seabird related global attributes found in DataFrame.header list
        
        Options
        -------
        
        Todo
        -----
        
        Retrieve PMEL header information from '@' comments in .cnv file or from separate
        header.txt file
        """        
        
        #From PMELheader
        
        self.rootgrpID.BAROMETER = Barometer
        self.rootgrpID.WIND_DIR = Wind_Dir
        self.rootgrpID.WIND_SPEED = Wind_Speed
        self.rootgrpID.AIR_TEMP = Air_Temp
        self.rootgrpID.WATER_DEPTH = Water_Depth
        self.rootgrpID.STATION_NAME = Station_Name
        self.rootgrpID.EPIC_FILE_GENERATOR = file_name
        self.rootgrpID.PROG_CMNT = Prog_Cmnt
        self.rootgrpID.EDIT_CMNT = Edit_Cmnt
        
        pass
        
    def dimension_init(self, len_time=1, len_depth=1, len_lat=1, len_lon=1 ):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time', 'depth', 'latitude', 'longitude']
        
        self.rootgrpID.createDimension( self.dim_vars[0], len_time ) #time
        self.rootgrpID.createDimension( self.dim_vars[1], len_depth ) #depth
        self.rootgrpID.createDimension( self.dim_vars[2], len_lat ) #lat
        self.rootgrpID.createDimension( self.dim_vars[3], len_lon ) #lon
        
        
    def variable_init(self):
        """
        U, V components of wind (m/s) - 'WU_422', 'WV_423'
        """

        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = ['T_25','ICEC_2025'], ['T','ICEC'], ["	SST (C)", "ICE CONC (%)"]
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = ['temp','icec'], ['',''], ['C','%'], [25,2025]
         


        rec_vars = ['time','time2','depth','latitude','longitude'] + rec_vars

        rec_var_name = ['', '', '', '', ''] + rec_var_name
        rec_var_longname = ['', '', '', '', ''] + rec_var_longname
        rec_var_generic_name = ['', '', '', '', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['f10.0', 'f10.0', 'f10.1', 'f10.4', 'f10.4'] + rec_var_FORTRAN
        rec_var_units = ['True Julian Day', 'msec since 0:00 GMT','m','degree_north','degree_west'] + rec_var_units
        rec_var_type= ['i4', 'i4'] + ['f4' for spot in rec_vars[2:]]
        rec_var_strtype= ['EVEN', 'EVEN', 'EVEN', 'EVEN', 'EVEN'] + ['' for spot in rec_vars[5:]]
        rec_epic_code = [624, 624,1,500,501] + rec_var_epic
        
        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[0]))#time2
        var_class.append(self.rootgrpID.createVariable(rec_vars[2], rec_var_type[2], self.dim_vars[1]))#depth
        var_class.append(self.rootgrpID.createVariable(rec_vars[3], rec_var_type[3], self.dim_vars[2]))#lat
        var_class.append(self.rootgrpID.createVariable(rec_vars[4], rec_var_type[4], self.dim_vars[3]))#lon
        
        for i, v in enumerate(rec_vars[5:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+5], rec_var_type[i+5], self.dim_vars))

        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, depth_level=-10., latitude=None, longitude=None, time1=None, time2=None, CastLog=False):
        """ """
        self.var_class[0][:] = time1
        self.var_class[1][:] = time2

        self.var_class[2][:] = depth_level #self.data[pressure_var].values
        self.var_class[3][:] = latitude
        self.var_class[4][:] = longitude #PMEL standard direction W is +

    def add_data(self, parameter, value):
        """ """

        di = self.rec_vars.index(parameter)
        self.var_class[di][:] = value
            
        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.History = self.History + ' ' + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                    + ' ' + new_history + '\n'
                    
    def close(self):
        self.rootgrpID.close()
    
class EPIC_NC_SST_cf(object):
    """ Class instance to generate a NetCDF file.  
    

    Standards
    ---------
    EPICNetCDF (PMEL) Standards  


    Usage
    -----
    
    Order of routines matters and no error checking currently exists
    ToDo: Error Checking
    
    Use this to create a nc file with all default values
        ncinstance = EPIC_NC()
        ncinstance.file_create()
        ncinstance.sbeglobal_atts()
        ncinstance.PMELglobal_atts()
        ncinstance.dimension_init()
        ncinstance.variable_init()
        ncinstance.add_coord_data()
        ncinstance.add_data()
        ncinstance.close()
    """ 
    
    
    nc_format = 'NETCDF3_CLASSIC'
    nc_read   = 'w'
    def __init__(self, savefile='data/test.nc', data=None):
        
        self.savefile = savefile
    
    def file_create(self):
        """ sets up filename, read/write permissions, netcdf format """
        
        rootgrpID = Dataset(self.savefile, EPIC_NC.nc_read, format=EPIC_NC.nc_format)
        self.rootgrpID = rootgrpID
        return ( rootgrpID )
        
    def sbeglobal_atts(self, DATA_CMNT='', coord_system="GEOGRAPHICAL", Water_Mass="G"):
        """
        Assumptions
        -----------
        
        Follows Mooring Format
        
        seabird related global attributes found in DataFrame.header list
        
        """
        
        self.rootgrpID.CREATION_DATE = datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")
        self.rootgrpID.DATA_CMNT = DATA_CMNT
        self.rootgrpID.COORD_SYSTEM = coord_system
        self.rootgrpID.WATER_MASS = Water_Mass
        
    def PMELglobal_atts(self, Barometer=9999, Wind_Dir=999, Wind_Speed=99,
                        Air_Temp=99.9, Water_Depth=9999, file_name='', Prog_Cmnt='', Edit_Cmnt='', Station_Name=''):
        """
        Assumptions
        -----------
        
        Format of DataFrame.name = 'dy1309l1_ctd001'
        
        seabird related global attributes found in DataFrame.header list
        
        Options
        -------
        
        Todo
        -----
        
        Retrieve PMEL header information from '@' comments in .cnv file or from separate
        header.txt file
        """        
        
        #From PMELheader
        
        self.rootgrpID.BAROMETER = Barometer
        self.rootgrpID.WIND_DIR = Wind_Dir
        self.rootgrpID.WIND_SPEED = Wind_Speed
        self.rootgrpID.AIR_TEMP = Air_Temp
        self.rootgrpID.WATER_DEPTH = Water_Depth
        self.rootgrpID.STATION_NAME = Station_Name
        self.rootgrpID.EPIC_FILE_GENERATOR = file_name
        self.rootgrpID.PROG_CMNT = Prog_Cmnt
        self.rootgrpID.EDIT_CMNT = Edit_Cmnt
        
        pass
        
    def dimension_init(self, len_time=1, len_depth=1, len_lat=1, len_lon=1 ):
        """
        Assumes
        -------
        Dimensions will be 'time', 'depth', 'lat', 'lon'
        
        Todo
        ----
        User defined dimensions
        """

        self.dim_vars = ['time', 'depth', 'latitude', 'longitude']
        
        self.rootgrpID.createDimension( self.dim_vars[0], len_time ) #time
        self.rootgrpID.createDimension( self.dim_vars[1], len_depth ) #depth
        self.rootgrpID.createDimension( self.dim_vars[2], len_lat ) #lat
        self.rootgrpID.createDimension( self.dim_vars[3], len_lon ) #lon
        
        
    def variable_init(self):
        """
        U, V components of wind (m/s) - 'WU_422', 'WV_423'
        """

        #build record variable attributes
        rec_vars, rec_var_name, rec_var_longname = ['T_25','ICEC_2025'], ['T','ICEC'], ["   SST (C)", "ICE CONC (%)"]
        rec_var_generic_name, rec_var_FORTRAN, rec_var_units, rec_var_epic = ['temp','icec'], ['',''], ['C','%'], [25,2025]
                


        rec_vars = ['time','depth','lat','lon'] + rec_vars

        rec_var_name = ['', '', '', ''] + rec_var_name
        rec_var_longname = ['', '', '', ''] + rec_var_longname
        rec_var_generic_name = ['', '', '', ''] + rec_var_generic_name
        rec_var_FORTRAN = ['', '', '', ''] + rec_var_FORTRAN
        rec_var_units = ['days since 1800-01-01 00:00:00','dbar','degree_north','degree_west'] + rec_var_units
        rec_var_type= ['f8'] + ['f4' for spot in rec_vars[1:]]
        rec_var_strtype= ['EVEN', 'EVEN', 'EVEN', 'EVEN'] + ['' for spot in rec_vars[4:]]
        rec_epic_code = [624,1,500,501] + rec_var_epic

        var_class = []
        var_class.append(self.rootgrpID.createVariable(rec_vars[0], rec_var_type[0], self.dim_vars[0]))#time1
        var_class.append(self.rootgrpID.createVariable(rec_vars[1], rec_var_type[1], self.dim_vars[1]))#depth
        var_class.append(self.rootgrpID.createVariable(rec_vars[2], rec_var_type[2], self.dim_vars[2]))#lat
        var_class.append(self.rootgrpID.createVariable(rec_vars[3], rec_var_type[3], self.dim_vars[3]))#lon
        
        for i, v in enumerate(rec_vars[4:]):  #1D coordinate variables
            var_class.append(self.rootgrpID.createVariable(rec_vars[i+4], rec_var_type[i+4], self.dim_vars))

        ### add variable attributes
        for i, v in enumerate(var_class): #4dimensional for all vars
            v.setncattr('name',rec_var_name[i])
            v.long_name = rec_var_longname[i]
            v.generic_name = rec_var_generic_name[i]
            v.FORTRAN_format = rec_var_FORTRAN[i]
            v.units = rec_var_units[i]
            v.type = rec_var_strtype[i]
            v.epic_code = rec_epic_code[i]
            
        self.var_class = var_class
        self.rec_vars = rec_vars

        
    def add_coord_data(self, depth_level=-10., latitude=None, longitude=None, time=None, CastLog=False):
        """ """
        self.var_class[0][:] = time

        self.var_class[1][:] = depth_level #self.data[pressure_var].values
        self.var_class[2][:] = latitude
        self.var_class[3][:] = longitude #PMEL standard direction W is +

    def add_data(self, parameter, value):
        """ """

        di = self.rec_vars.index(parameter)
        self.var_class[di][:] = value
            
        
    def add_history(self, new_history):
        """Adds timestamp (UTC time) and history to existing information"""
        self.History = self.History + ' ' + datetime.datetime.utcnow().strftime("%B %d, %Y %H:%M UTC")\
                    + ' ' + new_history + '\n'
                    
    def close(self):
        self.rootgrpID.close()
    


"""------------------------------------------------------------------------------------"""

def main():
    """ Nothing to do here """
    
if __name__ == "__main__":
    main() 