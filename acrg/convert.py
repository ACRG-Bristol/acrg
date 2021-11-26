import os
import json

from acrg.obs.read import synonyms
from acrg.config.paths import Paths

acrg_path = Paths.acrg

import numpy as np

def molar_mass(species):
    '''
    This function extracts the molar mass of a species from the acrg_species_info.json file.
    Returns:
        float : Molar mass of species
    '''
    species_info_file = os.path.join(acrg_path,"data/species_info.json")
    with open(species_info_file) as f:
            species_info=json.load(f)
    species_key = synonyms(species, species_info)
    molmass = float(species_info[species_key]['mol_mass'])
    return molmass

def mol2g(value,species):
    ''' Convert a value in moles to grams '''
    molmass = molar_mass(species)
    return value*molmass

def prefix(units):
        ''' Convert to a different prefix. Note there are probably libraries to do this. '''
        if units is None:
            unit_factor = 1.
        elif units == 'T':
            unit_factor=1.e12
        elif units == 'G': 
            unit_factor=1.e9
        elif units == 'P': 
            unit_factor=1.e15
        elif units == 'M': 
            unit_factor=1.e6
        else:
            print('Undefined prefix: outputting in g/yr')
            unit_factor=1.
        return unit_factor


def concentration(units):
    '''
    Conversion between mol/mol to parts-per- units
    '''
    unit_factor = 1e-12 if units.lower() == 'ppt' else \
                  1e-9 if units.lower()  == 'ppb' else \
                  1e-6 if units.lower()  == 'ppm' else \
                  1
    if unit_factor==1:
        print('Undefined prefix')
        
    return unit_factor

def convert_lons_0360(lons):
    '''
    Convert longitude values onto a 0-360 range from -180-180 range. Uses floored division. 
    
    Args:
        lons (arr):
            1D array of longitudes.            
    Returns:
        array:
            Longitudes on 0-360 range.           
    '''
 
    div = lons // 360

    return lons - div*360

def latlon_distance(lat, lon, units=None, verbose=True):
    '''
    Get the distance between 2 points in km
    
    Inputs
    ------
    lat, lon: list
        latitude and longitude of the 2 points, degrees
    units: str, optional
        units to return the distance in, defaults to km
        should be 'm', 'metres', 'meters', 'km', 'kilometres', 'kilometers'
    verbose: bool, optional
        default: True
    '''
    latlon = {key: [np.deg2rad(ll) for ll in latlon_deg]
              for key, latlon_deg in {'lat': lat, 'lon': lon}.items()}
    d_latlon = {key: abs(ll[1] - ll[0]) for key, ll in latlon.items()}

    radius_earth = 6378.137    # radius of the earth, km
    radius_earth = radius_earth * 1e3 if units in ['m', 'metres', 'meters'] else radius_earth

    # apply haversine formulae
    a = (np.sin(d_latlon['lat'] / 2))**2
    b = (np.sin(d_latlon['lon'] / 2))**2 * np.cos(latlon['lat'][0]) * np.cos(latlon['lat'][1])
    c = 2 * np.arcsinh(np.sqrt(a+b))
    return radius_earth * c

def convert_degrees_metres(lat, d_lat, d_lon, units=None, verbose=True):
    '''
    Convert the size of a pixel from degrees to km
    
    Inputs
    ------
    lat: list or int or float
        latitude of the centre of the pixel, degrees
    d_lat, d_lon: int or float
        size of the pixel, degrees
    units: str, optional
        units to return the distance in, defaults to km
        should be 'm', 'metres', 'meters', 'km', 'kilometres', 'kilometers'
    verbose: bool, optional
        default: True
    '''
    units = 'km' if units is None else units
    if verbose:
        print(f'Converting to units of {units}')
    
    if units not in ['m', 'metres', 'meters', 'km', 'kilometres', 'kilometers']:
        print('Please select units from: [m, metres, meters, km, kilometres, kilometers]')
        return None
    
    radius_earth = 6378.137    # radius of the earth, km
    radius_earth = radius_earth * 1e3 if units in ['m', 'metres', 'meters'] else radius_earth

    # convert all values to radians
    d_lat = np.deg2rad(d_lat)
    d_lon = np.deg2rad(d_lon)
    lat = np.deg2rad(lat)
 
    # get the 2 points of latitude
    lat = [lat - d_lat, lat + d_lat]
 
    # apply haversine formulae
    a = (np.sin(d_lat / 2))**2 + (np.sin(d_lon / 2))**2 * np.cos(lat[0]) * np.cos(lat[1])
    c = 2 * np.arcsinh(np.sqrt(a))
    return radius_earth * c