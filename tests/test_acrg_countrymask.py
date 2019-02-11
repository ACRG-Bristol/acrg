#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:13:44 2019

@author: rt17603
"""

import pytest
import acrg_countrymask as countrymask
import os

acrg_path = os.getenv("ACRG_PATH")

#%%
##### Create list of current NAME footprints ###################

@pytest.fixture(scope="module")
def name_domains():
    #domains = ['AUSTRALIA','CARIBBEAN','EASTASIA','EUROPE','NAMERICA','PACIFIC','SOUTHAFRICA','SOUTHASIA','WESTUSA']
    domains = ['AUSTRALIA','EASTASIA','EUROPE','NAMERICA']
    return domains

#%%
####### Tests for domain_volume function ######
def test_incorrect_domain():
    '''
    Test Exception is raised when incorrect domain is specified as a str
    '''
    with pytest.raises(Exception):# as e_info:
        domain = "MMM"
        countrymask.domain_volume(domain)

def test_wrong_input():
    '''
    Test Exception is raised when incorrect domain is specified as something other than a str
    '''
    with pytest.raises(Exception):# as e_info:
        domain = 1
        countrymask.domain_volume(domain)

def test_other_directory():
    '''
    Test that function can used with file read from a different directory
    '''
    fp_dir = os.path.join(acrg_path,"tests/files/NAME/fp/")
    print "acrg_path",acrg_path
    print "fp_dir",fp_dir
    domain = "EUROPE"
    fp_lat,fp_lon,fp_height = countrymask.domain_volume(domain,fp_directory=fp_dir)
    
    assert fp_lat is not None
    assert fp_lon is not None
    assert fp_height is not None

@pytest.mark.long
def test_all_domains(name_domains):
    '''
    Test all (current) domains are recognised and can be accessed
    '''
    for domain in name_domains:
        fp_lat,fp_lon,fp_height = countrymask.domain_volume(domain)

##########


#Samoa lat and lon 13.7590° S, 172.1046° W