#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:21:27 2020

@author: al18242

This file contains a method of obtaining the version of the code used to generate output.
"""

import subprocess

from acrg.config.paths import Paths
acrg_path = Paths.acrg

def code_version():
    '''   
    Use git describe to return the latest tag (and git hash if applicable).
    
    Returns
    -------
    version : String defining the version of the code used, or "Unknown" if git is unavailable

    '''
    
    try:
        output = subprocess.run(['git', 'describe'], capture_output=True,
                                cwd=acrg_path,universal_newlines=True)
        #remove newlines and cast as string
        version = str(output.stdout.strip('\n'))
    except:
        print("WARNING: Unable to identify version using git. Check that git is available to the python process.")
        #TODO: Add backup versioning method?
        version = "Unknown"
        
    return version