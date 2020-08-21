#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:21:27 2020

@author: al18242

This file contains a method of obtaining the version of the code used to generate output.
"""

import subprocess

from acrg_config.paths import paths
acrg_path = paths.acrg

def code_version():
    '''    
    Returns
    -------
    version : String defining the version of the code used

    '''
    
    try:
        output = subprocess.run(['git', 'describe'], capture_output=True,
                                cwd=acrg_path,universal_newlines=True)
        #remove newlines and cast as string
        version = str(output.stdout.strip('\n'))
    except:
        print("WARNING: Unable to identify version using git.")
        #TODO: Add backup versioning method?
        version = "Unknown"
        
    return version