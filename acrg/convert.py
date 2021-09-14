import os
import sys
import json
import acrg_obs

from acrg_config.paths import paths
acrg_path = paths.acrg

def molar_mass(species):
    '''
    This function extracts the molar mass of a species from the acrg_species_info.json file.
    Returns:
        float : Molar mass of species
    '''
    species_info_file = os.path.join(acrg_path,"acrg_species_info.json")
    with open(species_info_file) as f:
            species_info=json.load(f)
    species_key = acrg_obs.read.synonyms(species, species_info)
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