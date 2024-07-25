# -*- coding: utf-8 -*-
"""
Created on Thu May 2 16:10:38 2024

@author: bq24992

This module provides a set of functions for processing TCCON data.
"""

import glob
import os
import re
import random
import math
import numpy as np
import pandas as pd
import xarray as xray
import datetime as dt
from collections import OrderedDict
import itertools

import acrg.obs as acrg_obs
from .barometric import pressure_at_height
from acrg.countrymask import domain_volume
from acrg.config.paths import Paths
from acrg.satellite.gosat import gosat_output,add_history_attr,ds_check_internal_unique,extract_dates,define_pressure_levels


data_path = Paths.data

home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/TCCON/")
fp_directory = os.path.join(data_path,'LPDM/fp_NAME_test/')
obs_directory = os.path.join(data_path,'obs_test/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files_test") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure_test/")


def tccon_output_name(ds,site,max_level=17,use_name_pressure=False,pressure_domain=None,
                      pressure_base_dir=name_pressure_directory,
                      pressure_max_days=31,pressure_day_template=True,name_directory=name_csv_directory,
                      file_per_day=False,max_points=60,overwrite=False):
    '''
    '''
    # Note column mapping here is out:in values rather than in:out because some columns in output do not map to input data variables.
    col_mapping = OrderedDict([('ID_Level','lev'),
                               ('Time','time'),
                               ('x','longitude'),
                               ('dx','dlon'),
                               ('y','latitude'),
                               ('dy','dlat'),
                               ('z','pressure_levels'),
                               ('dz',None)])
    
    try:
        ds[col_mapping['dx']]
        ds[col_mapping['dy']]
    except KeyError:
        sounding = 10.5
        dlat = np.zeros(len(ds[col_mapping['y']]))+distance_lat(sounding)
        dlon = np.array([distance_lon(sounding,lat) for lat in ds[col_mapping['y']]])
        axis1 = col_mapping['Time']
        dlat_da = xray.DataArray(dlat,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
        dlon_da = xray.DataArray(dlon,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
        ds = ds.assign(**{col_mapping['dx']:dlon_da,col_mapping['dy']:dlat_da})
    
    columns=list(col_mapping.keys())
    
    if site is None:
        site = 'global'
    
    split_dim = "time"
    out_index = "ID_Level"
    
    name_pressure_convert = 100. # Input from GOSAT is in hPa and need to convert to Pa for NAME
    
    all_dates = extract_dates(ds,dim=split_dim)
    dates = np.unique(all_dates)
    
    pressure_levels,dpressure = define_pressure_levels(ds,pressure_domain=pressure_domain,
                                                       use_name_pressure=use_name_pressure,
                                                       pressure_base_dir=pressure_base_dir,
                                                       max_days=pressure_max_days,
                                                       day_template=pressure_day_template)
    
    network = "tccon"#find_network(site)[0] # Using first site as default
    
    name_directory = os.path.join(name_directory,network)
    if not os.path.isdir(name_directory):
        os.makedirs(name_directory)
    
    for date in dates:
        wh_date = np.where(all_dates == date)[0] # Find indices for each different date
        for ID,index in enumerate(wh_date):
            data = OrderedDict()
            ID_str = '{num:03d}'.format(num=ID+1) # Number for each point in a day
            
            for col in columns:
                name = col_mapping[col]
                if name:
                    ds_column = ds[name]
                    if name == "pressure_levels":
                        data[col] = pressure_levels[index][:max_level]*name_pressure_convert
                    elif name == "lev":
                        
                        if file_per_day:
                            # Write out ID as a str combination of the point number and level each with leading zeros
                            data[col] = np.array([ID_str+value.zfill(2) for value in (ds_column[:max_level].values+1).astype("str")])
                        else:
                            # Write out ID (levels) as two digit numbers with a leading zero e.g. 01, 02, ... 21, 22
                            data[col] = np.array([value.zfill(2) for value in (ds_column[:max_level].values+1).astype("str")])
                    else:
                        data[col] = np.array([ds_column[index].values]*max_level) # Populate each level with the same values e.g. lat,lon,dlat,dlon
                else:
                    if col == "dz":
                        data[col] = dpressure[index][:max_level]*name_pressure_convert

            df_output = pd.DataFrame(data,columns=columns)
            df_output.set_index(out_index,inplace=True)
            
            if file_per_day:
                filename = site+"_"+date.replace('-','')
                if max_points:
                    if len(wh_date) > max_points:
                        letter_split = chr(ord("A") + ID//max_points)
                        filename = "{}-{}".format(filename,letter_split)
                    ID_1 = ID%max_points
                else:
                    ID_1 = ID
                filename += '.csv'
                filename = os.path.join(name_directory,filename)
                
                if not overwrite:
                    if ID_1 == 0 and os.path.isfile(filename):
                        print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(filename))
                        break
                    elif not os.path.isfile(filename):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f',mode='a',header=False)    
                else:
                    if (not os.path.isfile(filename)) or (os.path.isfile(filename) and ID_1 == 0):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f',mode='a',header=False)
            else:
                filename = "{site}_{date}-{ID}.csv".format(site=site,date=date.replace('-',''),ID=ID_str)
                filename = os.path.join(name_directory,filename)
                if not overwrite:
                    if os.path.isfile(filename):
                        print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(filename))
                        break
                    else:
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                else:
                        print('Writing to filename: {}'.format(filename))    
                        df_output.to_csv(filename,float_format='%f')
                
    

def tccon_add_coords(ds,data_vars=[]):
########################## Begin change ##########################
    dims_current = ['time', 'prior_altitude', 'ak_altitude']
##########################  End change  ##########################
    dims = ['time','lev']
    
    if list(ds.dims.keys()) != dims_current and list(ds.dims.keys()) != dims_current[::-1]: # Check dimensions match what we expect (check forward and reverse list)
        print('WARNING: Do not recognise dimensions of input gosat dataset. Unable to add new dimensions.')
        return None
    #dims_current = ['n','m']
    
    # Define co-ordinate values. Extract from Dataset if already present (e.g. time), construct if not (e.g. lev)
    coords = {}
    for i,d in enumerate(dims):
        try:
            coords[d] = ds[d].values # Extract values for each coord from Dataset if present (e.g. time)
        except (AttributeError,KeyError):
            current_dim = dims_current[i]
            dim_length = ds.dims[current_dim]
            coords[d] = np.arange(dim_length) # If no values defined, create integer list based on dimension size
    
    if not data_vars:
        data_vars = ds.data_vars
    
    data = OrderedDict()
    for name in data_vars:
        if name not in dims:
            if len(ds[name].dims) == 1:
                if ds[name].values.size == coords[dims[0]].size:
                    data[name] = (dims[0],ds[name].values)
                elif ds[name].values.size == coords[dims[1]].size:
                    data[name] = (dims[1],ds[name].values)
                else:
                    raise ValueError("Cannot find a dimension of matching size.")
            elif len(ds[name].dims) == 2:
                data[name] = ([dims[0],dims[1]],ds[name].values)
########################## Begin change ##########################
    for name in dims_current:
        if name not in dims:
            if ds[name].values.size == coords[dims[0]].size:
                data[name] = (dims[0],ds[name].values)
            elif ds[name].values.size == coords[dims[1]].size:
                data[name] = (dims[1],ds[name].values)
##########################  End change  ##########################

    #coord_dict = {dims[0]:dim1,dims[1]:dim2}
    
    ds_new = xray.Dataset(data,coords=coords,attrs=ds.attrs)
    for dv in ds_new.data_vars:
        ds_new[dv].attrs = ds[dv].attrs
    for coord in ds_new.coords:
        try:
            ds_new[coord].attrs = ds[coord].attrs
        except (AttributeError,KeyError):
            if coord == "lev" or coord == "level":
                ds_new[coord].attrs["short_description"] = "Number for each level within the vertically resolved data."
                ds_new[coord].attrs["long_name"] = "level"
    
    # Add attribute describing modification made to original data
    mod_attr = "nc_restructure"
    ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Dataset rearranged to re-assign {} dimensions as {} coordinates, extracting values from dataset where present.".format(dims_current,dims).replace("'","")
    
    return ds_new
    
def tccon_process_file(filename,site,species="ch4",lat_bounds=[],lon_bounds=[],domain=None,
                       coord_bin=None,quality_filt=True,bad_pressure_filt=True,
                       name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],                       
                       mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                       pressure_domain=None,pressure_max_days=31.,pressure_day_template=True,
                       write_nc=False,output_directory=obs_directory,
                       write_name=False,name_directory=name_csv_directory,
                       file_per_day=False,max_name_level=17,max_name_points=None,overwrite=False,
                       verbose=True):
    '''
    '''
    print("\n\n\n############################################################")
    print("WARNING : Default max_name_level has to be changed (see what level is ~15km)")
    print("############################################################\n\n\n")  
    
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
    tccon = xray.open_dataset(filename)#[['time','prior_time',
                                        #  'lat','long','zobs','zmin',
                                        #  'prior_altitude',
                                        #  'ak_altitude','ak_pressure',
                                        #  'ak_xch4','prior_ch4','xch4','xch4_error',
                                        #  'airmass','tout','pout','hout',
                                        #  'prior_temperature','prior_pressure','prior_density']]
    
    tccon = tccon_add_coords(tccon)
    tccon = tccon.assign({'pressure_levels':(('time','lev'), \
                           tccon.ak_pressure.values[np.newaxis,:] \
                           *np.ones((tccon.time.values.size,tccon.lev.values.size)))
                            })
    
    print("\n\n\n############################################################")
    print("WARNING : Pressure_weights is crap (but I wanted to go on)")
    print("Should I put integrator operator in this variable?")
    tmp = tccon.ak_pressure.values[:-1]- tccon.ak_pressure.values[1:]
    tmp = np.concatenate([tmp,[tccon.ak_pressure.values[-1]]])/tccon.ak_pressure.values[0]
    tccon = tccon.assign({'pressure_weight':(('time','lev'), \
                        tmp*np.ones((tccon.time.values.size,tccon.lev.values.size)))
                        })
    print("############################################################\n\n\n")  

    print("\n\n\n############################################################")
    print("WARNING :  Equivalence of variables has to be checked") 
    print("('ak_xch4':'xch4_averaging_kernel','prior_ch4':'ch4_profile_apriori','xch4_error':'xch4_uncertainty')") 
    tccon = tccon.rename({'long':'longitude','lat':'latitude',
                          'ak_xch4':"xch4_averaging_kernel",
                          "prior_ch4":"ch4_profile_apriori",
                          "xch4_error":"xch4_uncertainty"})
    tccon = tccon.assign(exposure_id=lambda x: x.time.astype(int).astype(str))
    tccon = tccon.assign(retr_flag=lambda x: x.time.astype(bool))
    print("############################################################\n\n\n")  

    tccon = tccon.sortby(axis)
    tccon = ds_check_internal_unique(tccon,axis) # Check time values are unique and slightly modify if necessary
    
    # if quality_filt:
    #     if verbose:
    #         print('Applying filter based on quality flag')
    #     tccon = gosat_quality_filter(tccon)
    # if bad_pressure_filt:
    #     if verbose:
    #         print('Applying filter based on bad pressure flag')
    #     tccon = gosat_pressure_filter(tccon) # Removes points with any pressure values of -9999.99

    # if name_sp_filt:
    #     if len(tccon[axis].values) > 0:
    #         if not name_filters:
    #             raise Exception("If name_sp_filt=True, name_filters must be specified.")
    #         if verbose:
    #             print('Applying filters based on NAME surface pressure.')
            
    #         tccon = name_pressure_filter(tccon,filters=name_filters,pressure_domain=pressure_domain,
    #                                      cutoff=cutoff,layer_range=layer_range,pressure_base_dir=pressure_base_dir,
    #                                      max_days=pressure_max_days,day_template=pressure_day_template)
        
    columns=["latitude","longitude"]
    
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
        if verbose:
            print('Add virtual {0} degree bins to make it same format as satellite data'.format(coord_bin))
        tccon = tccon.assign({'dlat':(('time'), coord_bin[0]*np.ones(tccon.time.values.size)),
                              'dlon':(('time'), coord_bin[1]*np.ones(tccon.time.values.size))})
        # if verbose:
        #     print('Looking at area within latitude and longitude bounds: {0},{1}'.format(lat_bounds,lon_bounds))
        # tccon = latlon_filter(tccon,lat_bounds,lon_bounds,columns=columns)

    tccon.attrs["input_file"] = os.path.split(filename)[1]

    if len(tccon[axis].values) > 0:
        if write_nc:
            gosat_output(tccon,site=site,species=species,output_directory=output_directory,
                         file_per_day=file_per_day,overwrite=overwrite)
        
        if write_name:
            tccon_output_name(tccon,site=site,max_level=max_name_level,use_name_pressure=use_name_pressure,
                              pressure_domain=pressure_domain,pressure_base_dir=pressure_base_dir,
                              pressure_max_days=pressure_max_days,pressure_day_template=pressure_day_template,
                              name_directory=name_directory,file_per_day=file_per_day,
                              max_points=max_name_points,overwrite=overwrite)
    else:
        print('No points extracted from {} for specified parameters.'.format(filename))

    return tccon


def tccon_process(site,species="ch4",start=None,end=None,
                  input_directory=input_directory,
                  lat_bounds=[],lon_bounds=[],domain=None,coord_bin=None,
                  quality_filt=True,bad_pressure_filt=True,
                  name_sp_filt=False,name_filters=[],cutoff=5.,layer_range=[50.,500.],
                  mode=None,use_name_pressure=False,pressure_base_dir=name_pressure_directory,
                  pressure_domain=None,pressure_max_days=31,pressure_day_template=True,
                  write_nc=False,output_directory=obs_directory,
                  write_name=False,name_directory=name_csv_directory,file_per_day=False,
                  max_name_level=17,max_name_points=None,overwrite=False):
    '''
    '''

    
    if input_directory.find("$DATA_PATH"):
        input_directory = input_directory.replace("$DATA_PATH",str(data_path))
    
    search_str = os.path.join(input_directory,"*.nc")
    files = glob.glob(search_str)
    files.sort()

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
    
    if len(files) > 0:
        for i,filename in enumerate(files):
            if i == 0:
                verbose = True
            else:
                verbose = False

            ds = tccon_process_file(filename,
                                    site=site,
                                    species=species,
                                    lat_bounds=lat_bounds,
                                    lon_bounds=lon_bounds,
                                    coord_bin=coord_bin,
                                    quality_filt=quality_filt,
                                    bad_pressure_filt=bad_pressure_filt,
                                    max_name_level=max_name_level,
                                    name_sp_filt=name_sp_filt,
                                    name_filters=name_filters,
                                    cutoff=cutoff,
                                    layer_range=layer_range,
                                    mode=mode,
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
                                    verbose=verbose)
            
            if i == 0:
                tccon = ds
            else:
                tccon = gosat.merge(ds)
                
    else:
        print("No gosat files found within directory: {0} using search_str: {1} for dates {2} - {3}".format(input_directory,search_str,start,end))
        tccon = None
    
    return tccon       
    

                 
        
