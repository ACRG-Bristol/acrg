# -*- coding: utf-8 -*-
"""
Created on Thu May 2 16:10:38 2024

@author: bq24992

This module provides a set of functions for processing TCCON data.
"""

import glob
import os
import numpy as np
import pandas as pd
import xarray as xray
from collections import OrderedDict

from acrg.countrymask import domain_volume
from acrg.config.paths import Paths
from acrg.satellite.common import output_name,output,add_coords,ds_check_internal_unique

import logging

logger = logging.getLogger("acrg.satellite.tccon")

data_path = Paths.data

home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/TCCON/")
fp_directory = os.path.join(data_path,'LPDM/fp_NAME_test/')
obs_directory = os.path.join(data_path,'obs_test/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files_test") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure_test/")

siteList = ["BI, Bialystok, Poland",
            "BR, Bremen, Germany",
            "BU, Burgos, Philippines",
            "CI, California Institute of Technology, Pasadena, California, USA",
            "DB, Darwin, Australia",
            "DF, Armstrong Flight Research Center, Edwards, CA, USA",
            "ET, East Trout Lake, Canada",
            "EU, Eureka, Canada",
            "FC, Four Corners, NM, USA",
            "GM, Garmisch, Germany",
            "HF, Hefei, China",
            "HW, Harwell, UK",
            "IF, Indianapolis, Indiana, USA",
            "IZ, Izana, Tenerife, Spain",
            "JC, Jet Propulsion Laboratory 1, Pasadena, California, USA",
            "JF, Jet Propulsion Laboratory 2, Pasadena, California, USA",
            "JS, Saga, Japan",
            "KA, Karlsruhe, Germany",
            "LH, Lauder 1, New Zealand",
            "LL, Lauder 2, New Zealand",
            "LR, Lauder 3, New Zealand",
            "MA, Manaus, Brazil",
            "NI, Nicosia, Cyprus",
            "NY, Ny-Ålesund, Svalbard, Norway",
            "OC, Lamont, Oklahoma, USA",
            "OR, Orleans, France",
            "PA, Park Falls, Wisconsin, USA",
            "PR, LERMA, Sorbonne Unversite, Paris, France",
            "RA, Reunion Island, France",
            "RJ, Rikubetsu, Hokkaido, Japan",
            "SO, Sodankylä, Finland",
            "TK, Tsukuba, Ibaraki, Japan",
            "WG, Wollongong, Australia",
            "XH, Xianghe, China"]

def filter_and_resample(ds,var_to_keep,quality_filt):
    ds = ds[var_to_keep]
    if quality_filt:
        logger.info("Applying filter based on variable'extrapolation_flags_ak_xch4'.")
        ds = ds.where(abs(ds['extrapolation_flags_ak_xch4'])!=2)
    ds.dropna('time').sortby('time')
    tmp = ds.resample(time='h').mean(dim='time')
    tmp['xch4_maxError'] = ds['xch4_error'].resample(time='h').max(dim='time')
    del tmp['extrapolation_flags_ak_xch4']
    tmp = tmp.dropna('time')

    return tmp

def define_var_attrs(ds,method):
    ds['xch4_maxError'].attrs={'long_name':'xch4_maxError',
                               'description':'max of xch4_error on resampling period',
                               'unit':'ppm',
                               'vmin':str(ds.xch4_maxError.values.min()),
                               'vmax':str(ds.xch4_maxError.values.max())}
    ds['obs'].attrs={'long_name':'perturbed observation',
                     'unit':'ppm',
                     'vmin':str(ds.obs.values.min()),
                     'vmax':str(ds.obs.values.max())}
    ds['ak_footprint'].attrs={'vmin':str(ds.ak_footprint.values.min()),
                              'vmax':str(ds.ak_footprint.values.max())}
    if method=='integration_operator':
        ds['obs'].attrs['description']='xch4-prior_xch4-sum(ak_footprint*prior_ch4,dim="ak_altitude")'

        ds['ak_footprint'].attrs['long_name']='transformed ak using integration operator method'
        ds['ak_footprint'].attrs['description']='ak_xch4*integration_operator'
        ds.attrs['derivation_method'] = 'Integration operator'
        ds['dry_to_wet'].attrs={'long_name':'dry_to_wet',
                                'description':'ratio to use to convert from dry to wet mole fraction. Derived as 1/(1+prior_h2o)',
                                'unit':1,
                                'vmin':str(ds.dry_to_wet.values.min()),
                                'vmax':str(ds.dry_to_wet.values.max())}

    elif method=='pressure_weights':
        ds['obs'].attrs['description']='xch4-prior_xch4_from_dry-sum(ak_xch4*hj*prior_ch4,dim="ak_altitude")'

        ds['ak_footprint'].attrs['long_name']='ak derived using pressure weight method'
        ds['ak_footprint'].attrs['description']='see doc xxx'  
        ds.attrs['derivation_method'] = 'Pressure weight' 
        
        ds['wet_to_dry'].attrs={'long_name':'wet_to_dry',
                                'description':'ratio to use to convert from wet to wet dry fraction. Derived as 1/(1-prior_h2o)',
                                'unit':1,
                                'vmin':str(ds.wet_to_dry.values.min()),
                                'vmax':str(ds.wet_to_dry.values.max())}

    return ds

def tccon_process_file(filename,site,network,species="ch4",
                       start=None,end=None,lat_bounds=[],lon_bounds=[],domain=None,
                       coord_bin=None,method='pressure_weight',
                       quality_filt=True,
                       name_sp_filt=False,
                       use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31.,pressure_day_template=True,
                       write_nc=False,output_directory=obs_directory,
                       write_name=False,name_directory=name_csv_directory,
                       file_per_day=False,max_name_height=17,max_name_points=None,overwrite=False,
                       verbose=True):
    '''
    '''

    axis="time"
    
    if use_name_pressure or name_sp_filt:
        if domain and not pressure_domain:
            pressure_domain = domain
        elif not domain and not pressure_domain:
            raise Exception("To access NAME surface pressure files, pressure_domain (or domain) must be \
                            specified. Current pressure_domain={}.".format(pressure_domain))
    
    if lat_bounds:
        if len(lat_bounds) == 1 or len(lat_bounds) > 2:
            raise ValueError('Lat bounds must be specified as a two item iterable (e.g. list). Current value: {0}'.format(lat_bounds))
    
    if lon_bounds:
        if len(lon_bounds) == 1 or len(lon_bounds) > 2:
            raise ValueError('Lon bounds must be specified as a two item iterable (e.g. list). Current value: {0}'.format(lon_bounds))
    
    if coord_bin:
        if len(coord_bin) > 2:
            raise ValueError('Coordinate bin should be a one or two item iterable (e.g. list). Current value: {0}'.format(coord_bin))
    
    if verbose:
        print("========================")
        print('\nProcessing file: {0}\n'.format(filename))
    tccon = xray.open_dataset(filename)

    if start or end : 
        tccon = tccon.sel(time=slice(np.datetime64(start),
                                     np.datetime64(end)))
        if tccon.time.size==0:
            raise ValueError(f'No data found for {site} between {start} and {end}.')
    
    # Align units
    if tccon['prior_ch4'].units=='ppb' and tccon['prior_xch4'].units=='ppm' :
        tccon['prior_ch4'] = tccon['prior_ch4']*1e-3
        tccon['prior_ch4']['units']='ppm'

    var_to_keep = ['lat','long','obs','xch4','prior_xch4','prior_ch4','ak_footprint',
                   'extrapolation_flags_ak_xch4','xch4_error',
                   'ak_pressure','pressure_weight']
    
    if method=="integration_operator":
        logger.warning("Method 'pressure_weight' should be use to process TCCON data. "
                       "Bugs have been affecting variable 'integration_operator' used in this method. \n"
                       "Please be sure that you know what you do.")
        
        var_to_keep.append('dry_to_wet')
        # Select the levels
        tccon = tccon.sel(ak_altitude=slice(0,20)).sel(prior_altitude=slice(0,20))
        
        # Derive obs vector, ak for the footprints and dry to wet conv factor
        tccon['ak_footprint'] = tccon['ak_xch4'].interp(ak_altitude=tccon.prior_altitude)\
            *tccon['integration_operator']
        tccon['obs'] = tccon['xch4']-tccon['prior_xch4']\
            +(tccon['ak_footprint']*tccon['prior_ch4']).sum(dim='prior_altitude')
        tccon['dry_to_wet'] = 1/(1+tccon['prior_h2o'])

        # Define pressure_weights
        tccon['pressure_weight'] = tccon['integration_operator']
        tccon['pressure_weight'].attrs={'long_name':'pressure_weight',
                                       'description':'pressure weight (not sure that what is inside is what we expect, equal to integration_operator)',
                                       'unit':'1',
                                       'vmin':str(tccon.pressure_weight.values.min()),
                                       'vmax':str(tccon.pressure_weight.values.max())}

        # Filter the data and resample to hourly
        tccon = filter_and_resample(tccon,var_to_keep,quality_filt)
        # Define attributes
        tccon = define_var_attrs(tccon,method)

        logger.warning("Are footprint*flux_prior dry or wet mole fractions? Assumed as dry here.")
    
    elif method=="pressure_weight":
        
        var_to_keep.append('wet_to_dry')
        # Derive pressure thickness
        press = tccon.ak_pressure.values[:-1]- tccon.ak_pressure.values[1:]
        press = np.concatenate([press,[tccon.ak_pressure.values[-1]]])/tccon.ak_pressure.values[0]
        tccon = tccon.assign({'dpj':(('prior_altitude'),press)})

        # Derive pressure weight (hj), wet to dry conversion factor,
        # dry mole fraction of water (fdry_h2o) and prior dry xch4
        M_dryH2O,M_dryAir = 18.0153,28.9647
        tccon['wet_to_dry'] = 1/(1-tccon['prior_h2o'])
        tccon['fdry_h2o'] = tccon['prior_h2o']*tccon['wet_to_dry']
        tccon['hj'] = tccon['dpj']/(tccon['prior_gravity']
                                    *M_dryAir
                                    *(1+(tccon['fdry_h2o']
                                         *M_dryH2O/M_dryAir)))
        tccon['prior_xch4_from_dry'] = (tccon['hj']*tccon['wet_to_dry']*tccon['prior_ch4']
                                        ).sum(dim='prior_altitude') \
                                        /tccon['hj'].sum(dim='prior_altitude')
        
        # Select the levels
        tccon = tccon.sel(ak_altitude=slice(0,20)).sel(prior_altitude=slice(0,20))

        # Derive obs vector, ak for the footprints
        tccon['ak_footprint'] = tccon['ak_xch4'].interp(ak_altitude=tccon.prior_altitude)\
            *tccon['hj']/tccon['hj'].sum(dim='prior_altitude')
        tccon['obs'] = tccon['xch4']-tccon['prior_xch4_from_dry']\
            +(tccon['ak_footprint']*tccon['prior_ch4']).sum(dim='prior_altitude')

        # Define pressure_weights
        tccon['pressure_weight'] = tccon['hj']/tccon['hj'].sum(dim='prior_altitude')
        tccon['pressure_weight'].attrs={'long_name':'pressure_weight',
                                       'description':'pressure weight (not sure that what is inside is what we expect)',
                                       'unit':'1',
                                       'vmin':str(tccon.pressure_weight.values.min()),
                                       'vmax':str(tccon.pressure_weight.values.max())}
        
        # Filter the data and resample to hourly
        tccon = filter_and_resample(tccon,var_to_keep,quality_filt)
        # Define attributes
        tccon = define_var_attrs(tccon,method)
        
        logger.warning("Are footprint*flux_prior dry or wet mole fractions? Assumed as dry here.")
    
    tccon['ak_pressure'] = tccon['ak_pressure'].interp(ak_altitude=tccon.prior_altitude)
    del tccon['ak_altitude']
    tccon = add_coords(tccon,network='TCCON')
    
    tccon = tccon.rename({'long':'longitude','lat':'latitude',
                          'xch4':'xch4_tccon',
                          'obs':'xch4',
                          'ak_footprint':"xch4_averaging_kernel",
                          "prior_ch4":"ch4_profile_apriori",
                          "prior_xch4":"xch4_apriori",
                          "xch4_maxError":"xch4_uncertainty",
                          'ak_pressure':'pressure_levels'})
    
    tccon = tccon.assign(exposure_id=lambda x: x.time.astype(int).astype(str))
    tccon = tccon.assign(retr_flag=lambda x: x.time.astype(bool))

    tccon = tccon.sortby(axis)
    tccon = ds_check_internal_unique(tccon,axis) # Check time values are unique and slightly modify if necessary
            
    if domain and not (lat_bounds and lon_bounds):
        if verbose:
            print("Extracting latitude and longitude bounds from footprints associated with domain: {}".format(domain))
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if not lat_bounds:
        lat_bounds = [min(tccon.latitude.values),max(tccon.latitude.values)]
    if not lon_bounds:
        lon_bounds = [min(tccon.longitude.values),max(tccon.longitude.values)]
    
    if coord_bin:
        logger.info(f'Add virtual {coord_bin} degree bins to make it same format as satellite data')
        tccon = tccon.assign({'dlat':(('time'), coord_bin[0]*np.ones(tccon.time.values.size)),
                              'dlon':(('time'), coord_bin[1]*np.ones(tccon.time.values.size))})

    tccon.attrs["input_file"] = os.path.split(filename)[1]

    if len(tccon[axis].values) > 0:
        if write_nc:
            output(tccon,site=site,network=network,species=species,output_directory=output_directory,
                   file_per_day=file_per_day,overwrite=overwrite)
        
        if write_name:
            logger.warning("Units inconsistency between max_level (level) and max_name_height (km). TO CHECK")
            output_name(tccon,network=network,site=site,max_level=max_name_height,use_name_pressure=use_name_pressure,
                              pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                              pressure_max_days=pressure_max_days,pressure_day_template=pressure_day_template,
                              name_directory=name_directory,file_per_day=file_per_day,
                              max_points=max_name_points,overwrite=overwrite)
    else:
        print('No points extracted from {} for specified parameters.'.format(filename))

    return tccon


def tccon_process(site,network='TCCON',species="ch4",start=None,end=None,
                  input_directory=input_directory,
                  lat_bounds=[],lon_bounds=[],domain=None,coord_bin=None,
                  quality_filt=True,
                  name_sp_filt=False,
                  use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                  pressure_domain=None,pressure_max_days=31,pressure_day_template=True,
                  write_nc=False,output_directory=obs_directory,
                  write_name=False,name_directory=name_csv_directory,file_per_day=False,
                  max_name_height=17,max_name_points=None,overwrite=False):
    '''
    '''

    if network.upper()!='TCCON':
        raise ValueError('This function should be used only for TCCON data.')
    network = network.upper()
    
    if input_directory.find("$DATA_PATH"):
        input_directory = input_directory.replace("$DATA_PATH",str(data_path))
    
    if len(site)==2:
        site_long_name = np.array(siteList)[[x.split(',')[0].lower()==site.lower() for x in siteList]][0]
        site_short_name = site_long_name.split(',')[0].upper()
    else : 
        idx = [(site.lower() in x.lower()) for x in siteList]
        if sum(idx)!=1:
            raise ValueError(f'Please check site name. {sum(idx)} site(s) found {np.array(siteList)[idx]}.')
        site_long_name = np.array(siteList)[idx][0]
        site_short_name = site_long_name.split(',')[0].upper()
            
    
    search_str = os.path.join(input_directory,f"{site_short_name.lower()}*.nc")
    files = glob.glob(search_str)
    if len(files)!=1:
        raise ValueError(f'{len(files)} files found for site {site} ({site_long_name}). \n'
                         f'Please check input directory {input_directory} and be sure that data are downloaded from https://tccondata.org/ .')
    filename = files[0]

    if domain and not (lat_bounds and lon_bounds):
        lat,lon,height = domain_volume(domain)
        lat_bounds = [np.min(lat),np.max(lat)]
        lon_bounds = [np.min(lon),np.max(lon)]
    
    if (use_name_pressure or name_sp_filt) and domain and not pressure_domain:
        pressure_domain = domain
    
    if not write_nc and not write_name:
        print("\n**** WARNING: TCCON_PROCESS IS SET TO NOT OUTPUT ANYTHING TO FILE. ****")
        print("Please cancel this run and restart with write_nc or write_name ")
        print("parameters set to True if you want to write the processed output to disk.")
        print("*************************************************************************")
    
    if site!=site_short_name:
        site = site_short_name
        logger.info(f'Name of input site changed to "{site}".')

    ds = tccon_process_file(filename,
                            site=site,
                            network=network,
                            start=start,end=end,
                            species=species,
                            lat_bounds=lat_bounds,
                            lon_bounds=lon_bounds,
                            coord_bin=coord_bin,
                            quality_filt=quality_filt,
                            max_name_height=max_name_height,
                            name_sp_filt=name_sp_filt,
                            use_name_pressure=use_name_pressure,
                            pressure_base_dir=pressure_base_dir,
                            pressure_domain=pressure_domain,
                            pressure_max_days=pressure_max_days,
                            pressure_day_template=pressure_day_template,
                            write_nc=write_nc,
                            output_directory=output_directory,
                            write_name=write_name,
                            name_directory=name_directory,
                            file_per_day=file_per_day,
                            max_name_points=max_name_points,
                            overwrite=overwrite,
                            verbose=True)
    
    return ds       
    
