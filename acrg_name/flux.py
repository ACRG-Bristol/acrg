# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:13:48 2015

@author: chxmr
"""
from __future__ import print_function

from builtins import input
from builtins import str
from builtins import object
import numpy as np
import datetime as dt
import os
import netCDF4 as nc
import getpass
import acrg_time
import re
import numpy as np 
import xarray as xray
import collections as c
import sys
import pandas as pd

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    data_path = os.getenv("DATA_PATH") 
else:
    from acrg_config.paths import paths
    data_path = paths.data

output_directory = os.path.join(data_path,'LPDM/emissions/')

def write(lat, lon, time, flux, species, domain,
          source, title, prior_info_dict,
          regridder_used = 'acrg_grid.regrid.regrid_3D',
          copy_from_year = None, climatology = False, flux_comments = None,
          output_directory = output_directory):
    '''Write a flux file for emissions
    
    Args:
        lat (arr): 
            1D array of latitudes
        lon (arr): 
            1D array of longitudes
        time (numpy.datetime64): 
            Either an array (format from xarray) or a datetime index (format from pandas).
            See 'Creating datetime64 data' : http://xarray.pydata.org/en/stable/time-series.html
            If making a 'climatology' then this can be passed any dates but must be of same length 
            as time dimension of flux variable.
        flux (array):
            2D array of size [lat x lon x time]. Should be in mol/m2/s.
        species (str): 
            Species of interest.
        domain (str): 
            String of domain area
        source (str): 
            Sources in file. E.g. 'ff' (fossil fuel), 'agriculture'.
            source = None if the file contains all sources of a species. 
            If multiple sources: -source- is a chain of sources: 'waste-and-agriculture' (with hyphens between words).
        title (str): 
            Gives more information about what is in the file, e.g. Fossil Fuel CO2.
        prior_info_dict (dict): 
            {'NAME_OF_PRIOR' : ['VERSION','RAW RESOLUTION', 'REFERENCE']}
        regridder_used (str, optional): 
            regrid function used. Default is 'acrg_grid.regrid.regrid_3D'.
        copy_from_year (str, optional): 
            If the data is the same as another year but with a different timestamp give the original year here as a string
            Default is None
        climatology (bool, optional): 
            If the data is a climatology set this to True and give detail using flux_comments.
            Default is False
        flux_comments (str, optional): 
            Extra comments. Default is None.
        output_directory (str, optional): 
            Output directory. Default is 'data_path + LPDM/emissions/'.
    
    Returns:
        None
        Writes flux file to netcdf
        Naming convention of output file:
            The output filename depends on both whether a source is specified and if the output is
            for climatology.
            "species"-"source"_"domain"_"year".nc      (source specified, not climatology)
            "species"-climatology-"source"_"domain".nc (source specified, climatology)
            "species"-total_"domain"_"year".nc         (no source specified, not climatology) 
            "species"-climatology-total_"domain".nc    (no source specified, climatology)
            e.g. co2-ff_EUROPE_2012.nc
            e.g. co2-climatology-ff_EUROPE.nc
            e.g. ch4-total_SOUTHAMERICA_2010.nc
            e.g. ch4-climatology-total_SOUTHAMERICA.nc
        
    Example:
        flux.write(lat, lon, time, flux, 'ch4', 'EUROPE', '2012',
          comments = comments, title="2012_CH4_emissions_EUROPE") 
    
    Todo: 
        Add some error checking (e.g. check that domain is correct)
    '''
    
    print("WARNING: Make sure time stamp is start of time period (i.e. 1st of month\
            for monthly data or 1st January for yearly data).")
    print("WARNING: Make sure coordinates are centre of the gridbox.")
    print("WARNING: Make sure fluxes are in mol/m2/s.")

    if climatology == True:
        name_climatology = "-climatology"
    else:
        name_climatology = ""

    if source == None:
        file_source = species + name_climatology+ '-total'
        source_name = file_source
    else:
        file_source = species + name_climatology + '-' + source
        source_name = file_source
    
    file_source = file_source.lower()
    species = species.lower()  
        
    # Check that the flux is in the correct shape
    if np.shape(flux) != tuple((np.shape(lat)[0], np.shape(lon)[0], np.shape(time)[0])):
        print("Flux doesn't have dimensions lat x lon x time")
        print("Reshape your flux array and try again")
        return
        
    #Set climatology to year 1900
    if climatology == True:
        if len(time) == 1:
            time = [np.datetime64('1900-01-01')]
        elif len(time) == 12:
            time = np.arange('1900-01', '1901-01', dtype='datetime64[M]')
        else:
            sys.exit('Expecting either yearly or monthly climatology. Make sure time dimension is of size 1 or 12.')
    
    if type(time[0]) == np.datetime64:
        #time=time
        if isinstance(time,np.ndarray):
            time = time.astype(dtype="datetime64[ns]")
        elif isinstance(time,list):
            time = [t.astype("datetime64[ns]") for t in time]
        else:
            time = time
    else:
        sys.exit('Time format not correct, needs to be a list of type numpy.datetime64. A DatetimeIndex will not work.\
                 To convert a DatetimeIndex to correct format: time = [np.datetime64(i) for i in DatetimeIndex]')

    if not os.path.exists(output_directory):
        raise IOError("Unable to write file to specified output directory. Does not exist: {}.".format(output_directory))        

    #Open netCDF file
    year = pd.DatetimeIndex([time[0]]).year[0]
    if copy_from_year != None:
        ncname = os.path.join(output_directory,domain,f"{file_source}_{domain}_{year}_copy-from-{copy_from_year}.nc")
    elif climatology == True:
        ncname = os.path.join(output_directory,domain,f"{file_source}_{domain}.nc")
    else:
        ncname = os.path.join(output_directory,domain,f"{file_source}_{domain}_{year}.nc")

    if os.path.isfile(ncname) == True:
        answer = input("You are about to overwrite an existing file, do you want to continue? Y/N ")
        if answer == 'N':
            sys.exit()
        elif answer == 'Y':
            pass
    
    flux_attrs = {"source" : source_name,
                  "units" : 'mol/m2/s',
                  "species" : species} 
    
    lat_attrs = {"long_name" : "latitude",
                 "units" : "degrees_north",
                 "notes" : "centre of cell"}
    
    lon_attrs = {"long_name" : "longitude",
                 "units" : "degrees_east",
                 "notes" : "centre of cell"}


    glob_attrs = c.OrderedDict([("title",title),
                                ("author" , getpass.getuser()),
                                ("date_created" , np.str(dt.datetime.today())),
                                ("number_of_prior_files_used" , len(list(prior_info_dict.keys())))])

    for i, key in enumerate(prior_info_dict.keys()):
        prior_number = i+1
        glob_attrs['prior_file_' + str(prior_number)] = key
        glob_attrs['prior_file_' + str(prior_number)+'_version'] = prior_info_dict[key][0]
        glob_attrs['prior_file_' + str(prior_number)+'_raw_resolution']=prior_info_dict[key][1]
        glob_attrs['prior_file_' + str(prior_number)+'_reference']=prior_info_dict[key][2]
    
    glob_attrs["regridder_used"]= regridder_used
    
    if flux_comments != None:
        glob_attrs['comments'] = flux_comments
        if copy_from_year != None:
            glob_attrs['comments'] = f"Fluxes copied from year {copy_from_year}. {flux_comments}" 
    
    elif copy_from_year != None:
        glob_attrs['comments'] = f"Fluxes copied from year {copy_from_year}."

    flux_ds = xray.Dataset({'flux':(['lat','lon','time'], flux, flux_attrs)},
                              coords = {'lat' : lat,
                                        'lon' : lon,
                                        'time' : time},
                              attrs = glob_attrs)
    
    flux_ds.lat.attrs = lat_attrs
    flux_ds.lon.attrs = lon_attrs
    flux_ds.time.attrs['notes'] = "Start of time period"
    
    if not os.path.exists(os.path.join(output_directory,domain)):
        print(f"Creating {domain} subdirectory in output directory: {output_directory}")
        os.makedirs(os.path.join(output_directory,domain))

    flux_ds.flux.encoding = {'zlib':True}                        
    flux_ds.to_netcdf(ncname, mode='w')    


class EDGARread(object):
    def __init__(self, filename_of_EDGAR_emissions):

        f = nc.Dataset(filename_of_EDGAR_emissions, 'r')
    
        #Get grid
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get flux of species
#        variables = f.variables.keys()
#        species = str(variables[2])
        variables = [str(i) for i in list(f.variables.keys())]
        for i in variables:
            while i not in ['lat','lon','time']:
                species = i
                if species is not None:
                    break
        flux = f.variables[species][:,:]
        units = f.variables[species].units

        f.close()
        
        #Get year and datetime date of Edgar file
        filename = os.path.split(filename_of_EDGAR_emissions)[1]
        match = re.compile(r'_\d{4}_')
        m = match.findall(filename)
        if len(m) == 1:
            year = m[0].strip('_')
            date = dt.datetime.strptime(year, '%Y')
        elif len(m) > 1:
            print("Can't find correct date.")
            year = None
            date = None
        
        species = species.split('_')[1]

        self.lon = lon
        self.lat = lat
        self.flux = flux
        self.species = species
        self.units = units
        self.year = year
        self.date = date
        
