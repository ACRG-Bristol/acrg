#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:50:50 2019

@author: lw13938
"""
from builtins import str
from acrg_name.flux import write
import xarray as xr
import os
import sys

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    data_path = os.getenv("DATA_PATH") 
else:
    from acrg_config.paths import paths
    data_path = paths.data

data_path = str(data_path)

def test_write(tmpdir):
    with xr.open_dataset(data_path+'/LPDM/emissions/EUROPE/n2o-ocean_EUROPE_2009.nc') as load:
        ds = load.load()
    del ds.attrs['author']
    del ds.attrs['date_created']
    lat = ds.lat.values
    lon = ds.lon.values
    time = ds.time.values
    flux = ds.flux.values
    species = 'n2o'
    domain = 'EUROPE'
    source= 'ocean'
    title = ds.title
    prior_info_dict = {ds.prior_file_1 : [ds.prior_file_1_version,ds.prior_file_1_raw_resolution, ds.prior_file_1_reference]}     
    write(lat, lon, time, flux, species, domain,
          source, title, prior_info_dict,
          output_directory = str(tmpdir))
    with xr.open_dataset(str(tmpdir)+'/EUROPE/n2o-ocean_EUROPE_2009.nc') as load:
        testds = load.load()
    del testds.attrs['author']
    del testds.attrs['date_created']
    assert testds.equals(ds)
    assert testds.attrs == ds.attrs
    
