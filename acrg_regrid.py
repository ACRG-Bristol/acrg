# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:09:23 2015

regrid_emissions(filename_of_EDGAR_emissions,filename_of_footprint): Input the filename and path
of the EDGAR emissions you want to regrid and the filename and path of a
footprint that defines the correct domain and grid spacing for the new grid.

The output will be a .nc file of the interpolated emissions on the new grid.

It will be saved in /data/shared/NAME/emissions/DOMAIN-NAME/SPECIES_DOMAIN-NAME_YEAR.nc

@author: ew14860
"""

import acrg_name as name
import acrg_agage as agage
import netCDF4 as nc
import numpy as np
import scipy.interpolate as sci
import os.path
import json
import getpass
import acrg_time
import re
import datetime as dt
import iris


#Read in EDGAR data
class EDGARread:
    def __init__(self, filename_of_EDGAR_emissions):

        f = nc.Dataset(filename_of_EDGAR_emissions, 'r')
    
        #Get grid
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get flux of species
#        variables = f.variables.keys()
#        species = str(variables[2])
        variables = [str(i) for i in f.variables.keys()]
        for i in variables:
            while i not in ['lat','lon','time']:
                species = i
                if species is not None:
                    break
        flux = f.variables[species][:,:]
        units = f.variables[species].units

        f.close()

        self.lon = lon
        self.lat = lat
        self.flux = flux
        self.species = species
        self.units = units

#Find the domain name from the footprint file        
def find_domain(filename_of_footprint):
    filename = os.path.split(filename_of_footprint)[1]
    splitfile = filename.split('_')
    if len(splitfile) == 3:
        domain = splitfile[1]
    elif len(splitfile) != 3:
        print "Can't find domain name, make sure file name is in the correct format: site_domain_year.nc."
    return domain

#Find the year the EDGAR emissions data is for    
def find_EDGAR_year(filename_of_EDGAR_emissions):
    filename = os.path.split(filename_of_EDGAR_emissions)[1]
    match = re.compile(r'_\d{4}_')
    m = match.findall(filename)
    if len(m) == 1:
        year = m[0].strip('_')
        date = dt.datetime.strptime(year, '%Y')
    elif len(m) > 1:
        print "Can't find correct date."
    return date, year
    
#Regridding EDGAR data
def regrid(filename_of_EDGAR_emissions,
           filename_of_footprint,
           input_species_name = None):

#   Read in EDGAR grid
    read_ed = EDGARread(filename_of_EDGAR_emissions)
    lon_ed = read_ed.lon
    lat_ed = read_ed.lat
    flux_ed = read_ed.flux
    units = read_ed.units
    
    if input_species_name is None:
        species = read_ed.species.split('_')[1]
    elif input_species_name is not None:
        species  = input_species_name
    
    if lon_ed[0] > 0:
        lon_ed = lon_ed-180
        flux_ed0 = flux_ed.copy()
        flux_ed = np.zeros(np.shape(flux_ed0))
        flux_ed[:,0:(len(lon_ed)/2)]=flux_ed0[:,(len(lon_ed)/2):]
        flux_ed[:,(len(lon_ed)/2):]=flux_ed0[:,0:(len(lon_ed)/2)]
    
    elif lon_ed[0] < 0:
        lon_ed = lon_ed
        flux_ed = flux_ed
        
    coord_ed = (lat_ed, lon_ed)
    
#   Interpolate the EDGAR grid
    interp1 = sci.RegularGridInterpolator(coord_ed, flux_ed, method = 'linear')

#   Read in footprint grid
    read_fp = name.footprints(filename_of_footprint)
    lon_fp = read_fp.lon.values
    lat_fp = read_fp.lat.values
    
    X, Y = np.meshgrid(lat_fp, lon_fp)
    
    coord_fp = (X,Y)
    
#   Interpolate EDGAR emissions on footprint grid
    flux = interp1.__call__(coord_fp, method = 'linear')

    return lat_fp, lon_fp, flux.T, species, units
  
def iris_regrid(filename_of_EDGAR_emissions, filename_of_footprint):
    emissions = iris.load_cubes(filename_of_EDGAR_emissions)
    footprint = iris.load_cubes(filename_of_footprint)
#    Linear scheme
#    regridded_emi = emissions.regrid(footprint, iris.analysis.Linear())
#    Area-weighted scheme
    scheme = iris.analysis.AreaWeighted(mdtol=0.5)
    regridded_emi = emissions.regrid(footprint, scheme)
    
    return regridded_emi

# Find the correct unit converter to give flux in mol/m2/s
def unit_converter(current_units):
    if type(current_units) is not str:
        current_units = str(current_units)
    conversion = {"kg m-2 s-1" : {"alt": ["kg/m2/s"], "function": kg2mol}}
    units = agage.synonyms(current_units, conversion)
    function = conversion[units]['function']
    
    return function

#Converts kg/m2/s to mol/m2/s    
def kg2mol(kgflux, species_in):
    acrg_path=os.path.split(os.path.realpath(__file__))
    with open(acrg_path[0] + "/acrg_species_info.json") as f:
        species_info=json.load(f)
        
    species = agage.synonyms(species_in, species_info)
    mol_mass = float(species_info[species]['mol_mass'])
    molflux = kgflux*(1/(mol_mass*10**-3))
    
    return molflux

def mol2kg(molflux, species_in):
    acrg_path=os.path.split(os.path.realpath(__file__))
    with open(acrg_path[0] + "/acrg_species_info.json") as f:
        species_info=json.load(f)
        
    species = agage.synonyms(species_in, species_info)
    mol_mass = float(species_info[species]['mol_mass'])
    kgflux = molflux*(mol_mass*10**-3)
    
    return kgflux

#Write the new grid to a .nc file
def write(lat, lon, flux, species, domain, year, EDGAR_filename = None, footprint_filename = None):
    
    if type(year) == dt.date:
        time = year
        year = str(year.year)
    elif type(year) == str:
        time = dt.datetime.strptime(year, '%Y')
        year = year
    
    species = species.lower()    
    
    #Open netCDF file
    ncname = '/data/shared/NAME/emissions/%s/%s_%s_%s.nc' %(domain, species, domain, year)

    if os.path.isfile(ncname) == True:
        answer = raw_input("You are about to overwrite an existing file, do you want to continue? Y/N")
        if answer == 'N':
            return
        elif answer == 'Y':
            f=nc.Dataset(ncname, 'w')
   
    elif os.path.isfile(ncname) == False:
        f=nc.Dataset(ncname, 'w')
        
    #Create dimensions
    f.createDimension('time', len([time]))
    f.createDimension('lat', len(lat))
    f.createDimension('lon', len(lon))
    
    #Header attributes
    f.title = "EDGAR emissions"
    f.author = getpass.getuser()
    f.date_created = np.str(dt.datetime.today())
    if EDGAR_filename is not None:
        f.EDGAR_file_used = EDGAR_filename
    if footprint_filename is not None:
        f.footprint_file_used = footprint_filename
    f.notes = "Created using a linear interpolator in the function scipy.interpolate.RegularGridInterpolator"

    #Time variable
    time_seconds, time_reference = acrg_time.convert.time2sec(time)
    nctime=f.createVariable('time', 'i', ('time',))
    nctime[:]= time_seconds
    nctime.long_name='time'
    nctime.standard_name='time'
    nctime.units='seconds since ' + np.str(time_reference)
    nctime.calendar='gregorian'

    #Longitude variable
    nclon=f.createVariable('lon', 'd', ('lon',))
    nclon[:]=lon
    nclon.long_name='longitude'
    nclon.standard_name='longitude'
    nclon.units='degrees_east'
    
    #Latitude variable
    nclat=f.createVariable('lat', 'd', ('lat',))
    nclat[:]=lat
    nclat.long_name='latitude'
    nclat.standard_name='latitude'
    nclat.units='degrees_north'
    
    #Mole fraction variable
    ncflux=f.createVariable('flux', \
        'd', ('lat', 'lon', 'time'))
    ncflux[:,:,:]=flux
    ncflux.units='mol/m2/s'

    f.close() 

#The function that does everything
def regrid_emissions(filename_of_EDGAR_emissions,
                     filename_of_footprint,
                     input_species_name = None,
                     output_species_name = None):
    
    lat, lon, flux, species, units = regrid(filename_of_EDGAR_emissions,
                                            filename_of_footprint,
                                            input_species_name)
    converter = unit_converter(units)
    molflux = converter(flux, species)
    domain = find_domain(filename_of_footprint)
    date, year = find_EDGAR_year(filename_of_EDGAR_emissions)
    
    if output_species_name is None:
        output_spec = species
    elif output_species_name is not None:
        output_spec = output_species_name
        
    write(lat, lon, molflux, output_spec, domain, year,
          EDGAR_filename = filename_of_EDGAR_emissions,
          footprint_filename = filename_of_footprint)



