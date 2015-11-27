# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 17:58:39 2015

@author: chxmr
"""

import numpy as np
import json
from os import getenv
from os.path import join
import acrg_agage as agage

acrg_path = getenv("ACRG_PATH")
with open(join(acrg_path, "acrg_species_info.json")) as f:
    species_info=json.load(f)
    
def kg2mol(kg, species_in):
    '''
    Convert kg to mol
    '''
    species = agage.synonyms(species_in, species_info)
    mol_mass = float(species_info[species]['mol_mass'])
    return kg*(1./(mol_mass*10**-3))

def mol2kg(mol, species_in):
    '''
    Convert mol to kg
    '''
    species = agage.synonyms(species_in, species_info)
    mol_mass = float(species_info[species]['mol_mass'])
    return mol*(mol_mass*10**-3)
