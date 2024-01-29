#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:59:47 2018

This module is for creating country mask xarray.Dataset objects and netCDF files for new domains 
(or re-creating previous domains).

This is based on the regionmask module. This uses the Natural Earth database with a scale of 1:50m.

Output file name is of the form: "country_DOMAIN.nc" e.g. country_SOUTHAMERICA.nc

The domain can be extracted from current footprint files or latitude and longitude arrays can be specified.
Any country mask files already present within the output directory WILL NOT be overwritten. Rename or move 
the previous file if you wish to create a new file using this method.

Function to use for creating the country file is:
 - create_country_mask(...)
 
New function added to create country files that include marine territories (Exclusive 
Economic Zones). Function to do this is :
 - create_country_mask_eez(...)

@author: rt17603
"""
import regionmask
import iso3166
import numpy as np
import xarray as xr
from collections import OrderedDict
import glob
import getpass
import os
import cartopy
import cartopy.io.shapereader as shapereader
from shapely.prepared import prep
from shapely.geometry import Polygon
import geopandas as gpd
import pandas as pd
from scipy import stats

from acrg.convert import convert_lons_0360
from acrg.config.paths import Paths

data_path = Paths.data

fp_directory = os.path.join(data_path,'LPDM/fp_NAME/')
country_directory = os.path.join(data_path,'LPDM/countries/')

def domain_volume(domain,fp_directory=fp_directory):
    '''
    The domain_volume function extracts the volume (lat, lon, height) within a domain from a related footprint file.
    
    Args:
        domain (str) : 
            Domain of interest (e.g. one of 'AUSTRALIA', 'CARIBBEAN','EASTASIA','EUROPE','NAMERICA','PACIFIC',
            'SOUTHAFRICA','SOUTHASIA','WESTUSA')
        fp_directory (str, optional) : 
            fp_directory can be specified if files are not in the default directory. 
            Must point to a directory which contains subfolders organized by domain.
        
    Returns:
        xarray.DataArray (3): 
            Latitude, longitude, height
    '''
    directory = os.path.join(fp_directory,domain)
    listoffiles = glob.glob(os.path.join(directory,"*.nc"))
    if listoffiles:
        filename = listoffiles[0]
        print('Using footprint file: {0} to extract domain'.format(filename))
        with xr.open_dataset(filename) as temp:
            fields_ds = temp.load()
        
        fp_lat = fields_ds["lat"].values
        fp_lon = fields_ds["lon"].values
        fp_height = fields_ds["height"].values
    
        return fp_lat,fp_lon,fp_height     
    else:
        raise Exception('Cannot extract volume for domain: {1}. No footprint file found within {0}'.format(directory,domain))
        #return None

def range_from_bounds(bounds,step,include_upper=False):
    '''
    Create an array covering a range from a set of bounds. Choose whether to include upper
    bound or not.
    
    Args:
        bounds (list) :
            Two item list containing the upper and lower bounds of the range
        step (int/float) :
            Step size for the range
        include_upper (bool, optional) :
            Whether to include the upper bound.
            Default = False
    Returns:
        np.array:
            Array created using np.arange of range between lower and upper bounds.
    '''
    if include_upper:
        return np.arange(bounds[0],bounds[1]+step,step)
    else:
        return np.arange(bounds[0],bounds[1],step)

def country_name(code,allow_non_standard=True,supress_print=False):
    '''
    Extract country name based on ISO3166 standard.
    Note: allows non-standard values based on Natural Earth and Marine Regions databases to be included if 
    allow_non_standard is set to True.
    
    Args:
        code (str) :
            Can accept alpha2, alpha3 or numeric code to extract the country name.
        allow_non_standard (bool, optional) :
            Allow additional search for non-standard country names based on extra values 
            specifically from the Natural Earth and Marine Regions databases.
            Default = True.
        supress_print (bool, optional) :
            Whether to supress printing of a warning message when unable to find a match.
            Default = False.
    Returns:
        str:
            Country name (based on ISO3166 standard)
    '''
    non_standards_ne = {"DRC":u"Democratic Republic of the Congo",
                     "Dem. Rep. Congo":u"Democratic Republic of the Congo",
                     "SW":u"Swaziland","eSwatini":u"Swaziland",
                     "FL":u"Liechtenstein","Liechtenstein":u"Liechtenstein",
                     "KO":u"Kosovo","Kosovo":u"Kosovo"}
    non_standards_mr = {"Ascension":u"Saint Helena, Ascension and Tristan da Cunha",
                        "ASC":u"Saint Helena, Ascension and Tristan da Cunha",
                        "Tristan da Cunha":u"Saint Helena, Ascension and Tristan da Cunha",
                        "TAA":u"Saint Helena, Ascension and Tristan da Cunha",
                        "Clipperton Island":u"France","CPT":u"France"}
    
    non_standards = non_standards_ne
    non_standards.update(non_standards_mr)
    
    try:
        return iso3166.countries.get(code).name
    except KeyError:
        try:
            if not supress_print:
                print("Unable to derive country name from input: {}".format(code))
                print("WARNING: Looking through non-standard list")        
            return non_standards[code]
        except KeyError:
            if not supress_print:
                print("Unable to derive country name from non-standard list".format(code))
            return None

def country_alpha2(code,allow_non_standard=True,supress_print=False):
    '''
    Extract country two-letter code based on ISO3166 standard.
    Note: allows non-standard values based on Natural Earth and Marine Regions databases
    to be included if allow_non_standard is set to True.
    
    Args:
        code (str) :
            Can accept alpha3, numeric code or country name.
        allow_non_standard (bool, optional) :
            Allow additional search for non-standard country names based on extra values 
            specifically from the Natural Earth and Marine Regions databases.
            Default = True.
        supress_print (bool, optional) :
            Whether to supress printing of a warning message when unable to find a match.
            Default = False.
    
    Returns:
        str:
            Two-letter country code (based on ISO3166 standard)
    '''
    non_standards_ne = {"DRC":u"CD","Dem. Rep. Congo":u"CD","SW":u"SZ","eSwatini":u"SZ",
                     "FL":u"LI","Liechtenstein":u"LI"}
    non_standards_mr = {"Ascension":u"SH","ASC":u"SH","Tristan da Cunha":u"SH","TAA":u"SH",
                        "Clipperton Island":u"FR","CPT":u"FR"}
   
    # ASC - Ascension Islands - SH designates Saint Helena, Ascension and Tristan da Cunha (British Overseas Territory)
    # TAA - Tristan de Cunha - SH designates Saint Helena, Ascension and Tristan da Cunha (British Overseas Territory)
    # CPT - Clipperton Island - Uninhabited island in Eastern Pacific, oversees minor territory of France, included within FR country code.
    
    non_standards = non_standards_ne
    non_standards.update(non_standards_mr)
    
    try:
        return iso3166.countries.get(code).alpha2
    except KeyError:
        if allow_non_standard:
            try:
                if not supress_print:
                    print("Unable to derive country alpha2 code from input: {}".format(code))
                    print("WARNING: Having to look through non-standards")
                return non_standards[code]
                
            except KeyError:
                if not supress_print:
                    print("Unable to derive country alpha2 code from non-standard list")
                return None
        else:
            if not supress_print:
                print("Unable to derive country alpha2 code from input: {}".format(code))
            return None

def accepted_extra_countries(code):
    '''
    Additional country names with no equivalent within the ISO3166 standard but included (at the moment)
    within the Natural Earth database.
    
    e.g. Vatican, N. Cyprus, Kosovo, Somaliland
    
    Args:
        code (str) :
            Country name to check against list of additional accepted country names.
    
    Returns:
        str/None:
            returns the country name if it matches the list, None otherwise.
    '''
    accepted = [u"Somaliland",u"Vatican",u"N. Cyprus",u"Kosovo"]
    
    if code in accepted:
        return code
    else:
        return None

def manage_country_descriptor(countries):
    '''
    Manages the country names containing a comma and assumes the string after the comma contains the
    descriptor of the country e.g. "Republic Of". In most cases it is assumed that this detail is not 
    needed in the final country name. Exceptions are detailed below.
    
    Rules for commas are as follows:
     - If the country requires the designation after a comma the name will be rearranged. At the moment,
     this includes "congo" and "korea".
     e.g. "Korea, Republic of" becomes "Republic of Korea"
     e.g. "Korea, Democratic People's Republic of" becomes "Democratic People's Republic of Korea"
     - Otherwise, any designation specified after a comma will be removed.
     e.g. "Venezuela, Bolivarian Republic of" becomes "Venezuela"
    
    WARNING: Changes countries list in place
    
    Args:
        countries (list) :
            List of country names to manage.
    
    Returns:
        list:
            Country names altered as described above where appropriate.
    '''
    descriptor_needed = ["congo","korea"]
    for i,name in enumerate(countries):
        if name.find(',') != -1:
            name_split = name.split(',')
            if len(name_split) == 2:
                country,descriptor = name_split
                if country.lower() in descriptor_needed:
                    countries[i] = descriptor.lstrip() + ' ' + country.rstrip()
                else:
                    countries[i] = country.rstrip()

    return countries
    

def region_mapping(upper=True):
    '''
    Extract mapping from the region number to the country name from the regionmask module.
    
    See manage_country_descriptor() function for how names containing commas are treated.
    
    Args:
        upper (bool, optional):
            Return country names all in upper case. Default = True.
    
    Returns:
        dict:
            Dictionary containing the mapping of the region number to the country name.
    '''
    all_regions = regionmask.defined_regions.natural_earth_v5_0_0.countries_50.region_ids
    #codes = [key for key in all_regions.keys() if isinstance(key,unicode) and len(key) == 2]
    #names = [country_name(code,supress_print=True) for code in codes]
    
    #names = [name.split(',')[0] if name else name for name in names]
    
    names = regionmask.defined_regions.natural_earth_v5_0_0.countries_50.names
    numbers = [all_regions[key] for key in names]

    # Check all numbers are covered
    max_num = np.max(numbers)
    for i in range(max_num+1):
        if i not in numbers:
            #for key,item in list(all_regions.items()):
            for key,item in all_regions.items():
                if item == i:
                    if isinstance(key,str) and len(key) > 1:
                        name = country_name(key,supress_print=True)
                        if name:
                            numbers.append(i)
                            names.append(name)
                            break
                        else:
                            if accepted_extra_countries(key):
                                numbers.append(i)
                                names.append(key)
                                break
    
    names = manage_country_descriptor(names)
    if upper:
        names = np.core.defchararray.upper(names)
    #region_dict = {num:name for num,name in list(zip(numbers,names))}
    region_dict = {num:name for num,name in zip(numbers,names)}

    return region_dict

def mask_index_reset(mask,regions=None):
    '''
    Reset index within mask from numbers derived from region_mask module to incrementally increasing values.
    
    Args:
	    mask (xarray.DataArray) :
	    	DataArray containing grid (nlat x nlon)
	    regions (list/None, optional):
	    	List of region numbers. If not specified, this will be extracted from mask.
	    	Default = None.
    
    Returns:
	    xarray.DataArray:
	    	Mask with values renumbered from 0-nregions
    '''
    
    if regions is None:
        regions = np.unique(mask)
    
    ## Normalise numbers in mask- have to add arbitrary number to avoid clashes when numbers are being reassigned
    add_num = 100000
    mask.values = mask.values+add_num
    for i,region_num in enumerate(regions):
        mask_num = region_num+add_num
        mask.values[np.where(mask == mask_num)] = i
    
    return mask

    
    
def create_country_mask(domain,lat=None,lon=None,reset_index=True,ocean_label=True,use_domain_extent=False,
                        fp_directory=fp_directory,write=True,output_dir=country_directory):
    '''
    Creates a country mask for the latitude and longitude range. Derives country data from
    Natural Earth database (at 1:50m resolution).
    
    If write=True, writes out file of the form "country_DOMAIN.nc" e.g. country_SOUTHAMERICA.nc
    Note: Will not overwrite a pre-existing file.
    
    Args:
        domain (str) :
            Domain associated with lat and lon grid.
            If lat and lon arrays are not specified, domain will be used to extract latitude and
            longitude grids.
        lat (np.array, optional) :
            Latitude grid values.
            If not specified, will attempt to extract latitude values from domain footprint file.
        lon (np.array, optional) :
            Longitude grid values.
            If not specified, will attempt to extract longitdue values from domain footprint file.
        reset_index (bool, optional) :
            Reset index values within grid to 0 - num regions rather than retain their
            original values within the Natural Earth input (regionmask input).
            May be useful to not reset the index when needing to compare to other inputs.
            Should be set to True when writing to file to match with other ACRG code for looking at 
            country totals.
            Default = True.
        ocean_label (bool, optional) :
            Add additional explicit ocean label as well as countries within domain.
            Default = True.
        use_domain_extent (bool, optional) :
            If lat and lon values are specified, still create region covering the domain area but
            mask out any values outside the lat, lon bounds (treat lat,lon as a sub-domain).
            Note: Only the lat and lon outer values will be used as the latitude and longnitude values
            in the grid must match the domain values.
        fp_directory (str, optional) :
            Base footprint directory to use to extract domain values. Uses footprint directory on data path
            by default and expects sub-directories of domain name. Only used if lat, lon are not specified.
        write (bool, optional) :
            Write produced dataset to file of the form "country_"DOMAIN".nc".
            Default = True.
        output_dir (str, optional) :
            Directory for writing output.
    
    Returns:
        xarray.Dataset:
            Dataset containing the country map
        
        If write:
            Writes the dataset to file.
    '''
    
    database = "Natural_Earth"
    scale = "1:50m"
    
    if lat is None and lon is None:
        lat,lon,height = domain_volume(domain, fp_directory=fp_directory)
    elif lat is None or lon is None:
        raise Exception("Latitude and Longitude arrays must both be specified. Otherwise domain can be used to find these values.")
    elif use_domain_extent:
        sub_lat = lat[:]
        sub_lon = lon[:]
        lat,lon,height = domain_volume(domain, fp_directory=fp_directory)
    
    if database == "Natural_Earth" and scale == "1:50m":
        # moves longitude values onto 0-360 range if lons have values less than 0 and greater than 180, as regionmask cannot use a lon range that has both 
        if any(lon < 0) & any(lon > 180):
            lon = convert_lons_0360(lon)
        mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_50.mask(lon,lat)
    elif database == "Natural_Earth" and scale == "1:110m":
        # moves longitude values onto 0-360 range if lons have values less than 0 and greater than 180, as regionmask cannot use a lon range that has both 
        if any(lon < 0) & any(lon > 180):
            lon = convert_lons_0360(lon)
        mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(lon,lat,xarray=True)

    if use_domain_extent:
        lat_outer = np.where((mask["lat"].values < sub_lat[0]) | (mask["lat"].values >= sub_lat[-1]))[0]
        lon_outer = np.where((mask["lon"].values < sub_lon[0]) | (mask["lon"].values >= sub_lon[-1]))[0]
        
        mask.values[lat_outer,:] = np.nan
        mask.values[:,lon_outer] = np.nan
        
        if ocean_label:
            print("Cannot create ocean label when using domain extent with a sub-domain")
            ocean_label = False
    
    ## Find region numbers and associate with countries    
    mask_flat = mask.values.flatten()
    regions = np.unique(mask_flat[~np.isnan(mask_flat)])
    region_dict = region_mapping()
    countries = [region_dict[region] for region in regions]

    ## Check for multiple region numbers with the same region name
    unique,inverse,counts = np.unique(countries,return_inverse=True,return_counts=True)
    if counts[counts>1].any():
        repeats = np.where(counts>1)[0]
        repeat_indices = [np.where(inverse == r)[0] for r in repeats]
        repeat_regions = [regions[ri] for ri in repeat_indices]
        
        for region_mult in repeat_regions:
            region_initial = region_mult[0]
            region_repeat = region_mult[1:]
            for region in region_repeat:
                mask.values[np.where(mask == region)] = region_initial
        
        index_remove = [ri for ri_mult in repeat_indices for ri in ri_mult[1:]]
        index_remove.reverse()
        for ri in index_remove:
            countries.pop(ri)
        regions = np.delete(regions,index_remove)

    ## Add ocean label if relevant and set appropriate ocean_ref
    if ocean_label:
        ocean_ref = -1
        countries.insert(0,"OCEAN")
        regions = np.insert(regions,0,[ocean_ref])
    else:
        ocean_ref = -(2)**31
    
    ## Set NaN values as ocean reference but make sure to keep any regions already labelled as 0.0 as 
    # nan_to_num function changes the nan values to 0.
    temp_num = 100000
    mask.values[np.where(mask == 0)] = temp_num
    mask.values = np.nan_to_num(mask)
    mask.values[np.where(mask == 0)] = ocean_ref
    mask.values[np.where(mask == temp_num)] = 0
    
    if reset_index:
        mask = mask_index_reset(mask,regions=regions)
   
    ## Turn mask output into a dataset and add additional parameters and attributes
    
    ds = mask.to_dataset(name="country")

    countries = np.array(countries).astype(str)
    if reset_index:
        ds["name"] = xr.DataArray(countries,dims="ncountries")
    else:
        ds["name"] = xr.DataArray(countries,coords={"ncountries":regions},dims={"ncountries":len(regions)})
    
    ds["country"].attrs = {"long_name":"Country indices"}
    ds["lat"].attrs = {"long_name":"latitude","units":"degrees_north"}
    ds["lon"].attrs = {"long_name":"longitude","units":"degrees_east"}
    ds["name"].attrs = {"long_name":"Country names"}
    
    attributes = OrderedDict([])
    attributes["title"] = "Grid of country extent across {} domain".format(domain)
    attributes["author"] = "File created by {}".format(getpass.getuser())
    attributes["database"] = "{} database with scale {} used to create this file. Python module regionmask (https://regionmask.readthedocs.io/en/stable/) used to extract data.".format(database,scale)
    attributes["extent"] = "Domain beween latitude {} - {}, longitude {} - {} with resolution {:.3f},{:.3f} degrees".format(np.min(lat),np.max(lat),np.min(lon),np.max(lon),lat[1]-lat[0],lon[1]-lon[0])

    ds.attrs = attributes
    
    if write:
        filename = "country_{}.nc".format(domain)
        filename = os.path.join(output_dir,filename)
        if not os.path.isfile(filename):    
            print('Writing output to {}'.format(filename))
            ds.to_netcdf(filename)
        else:
            print("ERROR: DID NOT WRITE TO FILE: {} already exists".format(filename))
    
    return ds

def mask_fill_gaps(mask_array,lat,lon):
    
    #TODO: create a test for fill_gaps
    print('\nFilling in gaps between masks by checking for grid cells where surrounding cells are indentical.')
    print('WARNING: This is experimental, when using this mask to estimate country emissions from countries '+
          'with complex coastlines you should manually check the mask before use.')
    print('WARNING: This does not work for any gaps that are on the edge of the domain.')
        
    # loop through grid to check for missing cells

    for i in np.arange(1,mask_array.shape[0]-1):
        for j in np.arange(1,mask_array.shape[1]-1):

            if mask_array[i,j] == 0.:
                    
                # find all surrounding grid cells
                surrounding = np.delete(mask_array[i-1:i+2,j-1:j+2].flatten(),4)

                # if all surrounding cells are identical, fill in the gap
                if np.all(surrounding == surrounding[0]) and surrounding[0] != 0.:
                    print('Filling gap at '+ str(lat[i]) + ' ' + str(lon[j]))

                    mask_array[i,j] = surrounding[0]

                # if all but one of the surrounding cells are identical, fill in the gap
                # (fills in gaps of two adjacent grid cells)
                elif np.count_nonzero(surrounding == stats.mode(surrounding)[0]) == 7. and stats.mode(surrounding)[0] != 0.:
                    print('Filling gap at '+ str(lat[i]) + ' ' + str(lon[j]))
                    mask_array[i,j] = stats.mode(surrounding)[0]
                
                # fills gaps where 6 of the surrounding celss are identical
                elif np.count_nonzero(surrounding == stats.mode(surrounding)[0]) == 6. and stats.mode(surrounding)[0] != 0.:
                    print('Filling gap at '+ str(lat[i]) + ' ' + str(lon[j]))
                    mask_array[i,j] = stats.mode(surrounding)[0]
                    
    return mask_array

def create_country_mask_eez_v12(domain,lat=None,lon=None,include_land_territories=True,
                            include_ocean_territories=True,reset_index=True,fill_gaps=True,
                            fp_directory=fp_directory,lat_lon_mask=False,sub_lats=None,sub_lons=None,
                            include_countries=None,include_uk_dn=False,
                            output_path=None,save=False):
    """
    
    NOTE: THIS FUNCTION IS CALLED CREATE_COUNTRY_MASK_EEZ_V12 IN THE ACRG REPO
    
    Creates a mask for all countries within the domain (or lat/lon bounds).
    Uses Natural Earth 10m land datasets and Admin_0_map_units datasets
    for specifiying country and ocean boundaries.
    Option to include EEZ (Exclusive Economic Zones e.g. marine 
    territories) is True by default. Remaining areas are given a value of 0.
    
    Based on code written by Daniel Hoare and Rachel Tunnicliffe.
    
    Args:
        domain (str):
            The domain the mask should cover, e.g. 'EUROPE'. If lat and lon are not specified,
            will use the domain to extract the latitudes and longitudes.
        lat (numpy array) (optional):
            Latitudes to create mask over. If not specified, will extract lats from a domain.
        lon (numpy array) (optional):
            Longitudes to create mask over. If not specified, will extract lons from a domain.
        include_land_territories (bool):
            If True, include land territories for all countries in the 
            country mask.
        include_ocean_territories (bool):
            If True, include areas within EEZ (marine territories) in the 
            country mask.
        reset_index (bool):
            If True, start indexing countries from 1 (with ocean as 0). 
            If False, uses the Natural Earth index values (as extracted by geopandas).
        fill_gaps (bool):
            Fills in gaps between the land and ocean masks that are missed and some 
        fp_directory (str, optional) :
            Base footprint directory to use to extract domain values. Uses footprint directory on data path
            by default and expects sub-directories of domain name.
        lat_lon_mask (bool):
            If True, masks all cells outside sub_lats and sub_lons (gives them values of np.nan).
        sub_lats (list or array) (optional):
            Either a list of [min_lat,max_lat] or a numpy array of lats. Any cells beyond these limits are masked.
        sub_lons (list or array) (optional):
            Either a list of [min_lon,max_lon] or a numpy array of lons. Any cells beyond these limits are masked.
         out_path (str) (optional):
            Path and name to save mask to. If None, does not save dataset.
    Returns:
        ds (xarray dataset):
            Country mask dataset.
    """
    
    print("If you are having issues with downloading the required files from natural earth:")
    print("Try installing cartopy version 0.20.0, this may fix the issue.")
   
    if lat is not None:
        lat_range = [lat[0],lat[-1]]
        lon_range = [lon[0],lon[-1]]
        print(f'\nCreating mask between lats: {lat_range} and lons: {lon_range}.')
        lats = lat
        lons = lon
    else:
        print(f'\nExtracting lats and lons from domain: {domain}.')
        # extract lats and lons from fp file
        lats,lons,heights = domain_volume(domain,fp_directory=fp_directory)

    # load in land files
    shpfilename_land = shapereader.natural_earth('10m','cultural','admin_0_map_units')

    df_land = gpd.read_file(shpfilename_land)
    
    print('Creating land mask...')
    mask = regionmask.mask_3D_geopandas(df_land['geometry'],lons,lats)
    
    land_regions = mask.region.values
    
    # load in ocean files
    if include_ocean_territories:

        shpfilename_ocean = os.path.join(data_path,'World_shape_databases/MarineRegions/eez_v12.shp')

        df_ocean = gpd.read_file(shpfilename_ocean)
        
        print('Creating ocean mask...')
        mask_ocean = regionmask.mask_3D_geopandas(df_ocean['geometry'],lons,lats)
        
        ocean_regions = mask_ocean.region.values
        
    else:
        print('Not including ocean territories')
        
    all_grid = np.zeros((lats.shape[0],lons.shape[0]))
    country_names = []
    country_codes = []
    
    if reset_index == True:
    
        #land_codes_all = np.unique(df_land.loc[land_regions]['ADM0_A3'])
        land_codes_all = pd.unique(df_land.loc[land_regions]['ADM0_A3'])
        
        land_codes_all = np.sort(land_codes_all)
        
    else:
    
        land_codes_all = []
        indexes_all = []
        
        for i,l_code in enumerate(df_land.loc[land_regions]['ADM0_A3']):
            if l_code not in land_codes_all:
                land_codes_all.append(l_code)
                indexes_all.append(df_land.loc[land_regions]['ADM0_A3'].index[i])

    if type(include_countries) == list or type(include_countries) == np.ndarray:
        print(f'Only including: {include_countries}.')
        land_codes_all = include_countries.copy()
        if reset_index == False:
            indexes_all = []
            for c in land_codes_all:
                indexes_all.append(np.where(df_land.loc[land_regions]['ADM0_A3'] == c)[0][0])
    else:
        print('Including all countries in land and ocean masks.')
            
    count = 1

    # loop through all countries
    for i,c in enumerate(land_codes_all):
        if reset_index == True:
            index_value = count
        else:
            print('WARNING - not resetting index will mean that the isle of man is not correctly included in the UK mask')
            index_value = indexes_all[i]

        country_names.append(df_land[df_land['ADM0_A3'] == c]['ADMIN'].values[0].upper())
        country_codes.append(c)
        
        if include_land_territories:
        
            # find shapefile dataset indexes that correspond to the country code 
            land_region_indexes = np.where(df_land['ADM0_A3'] == c)[0]

            for r in land_region_indexes:

                if r in land_regions:

                    # find mask index that corresponds to the dataset index
                    land_mask_index = np.where(land_regions == r)[0][0]

                    all_grid[np.where(mask[land_mask_index].values == True)] = index_value
                    
        elif i == 0:
            
            print('Not including land territories')

        if include_ocean_territories:    
            
            # find shapefile dataset indexes that correspond to the country code
            ocean_region_indexes = np.where(df_ocean['ISO_TER1'] == c)[0]

            for r in ocean_region_indexes:
                
                # checks that ocean region is in domain and not disputed territory
                if r in ocean_regions:
                    try:
                        np.isnan(df_ocean['ISO_TER2'][r])
                    
                        # find ocean mask index that corresponds to the dataset index
                        ocean_mask_index = np.where(ocean_regions == r)[0][0]

                        all_grid[np.where(mask_ocean[ocean_mask_index].values == True)] = index_value
                    except:
                        print(f'Ocean region {r} is disputed territory.')

        count += 1
        
    # fill UK crown dependencies with the UK value
    if 'GBR' in country_codes and 'IMN' in country_codes:
        
        print('\nAdding the Isle of Man (IMN), Guernsey (GGY) and Jersey (JEY) to the UK mask\n')
        
        imn_index = np.where(all_grid == np.where(np.array(country_codes) == 'IMN')[0][0]+1)
        gbr_value = np.where(np.array(country_codes) == 'GBR')[0][0]+1
        
        all_grid[imn_index] = gbr_value
        
        
    if 'GBR' in country_codes and 'JEY' in country_codes:
        
        jey_index = np.where(all_grid == np.where(np.array(country_codes) == 'JEY')[0][0]+1)
        gbr_value = np.where(np.array(country_codes) == 'GBR')[0][0]+1

        all_grid[jey_index] = gbr_value
        
    if 'GBR' in country_codes and 'GGY' in country_codes:
        
        ggy_index = np.where(all_grid == np.where(np.array(country_codes) == 'GGY')[0][0]+1)
        gbr_value = np.where(np.array(country_codes) == 'GBR')[0][0]+1

        all_grid[ggy_index] = gbr_value
    
    if fill_gaps:
        
        all_grid = mask_fill_gaps(all_grid,lats,lons)
        all_grid = mask_fill_gaps(all_grid,lats,lons)
        
    
    if lat_lon_mask == True:
        print(f'Masking all values beyond lats: {sub_lats[0]}:{sub_lats[-1]} and lons: {sub_lons[0]}:{sub_lons[-1]}')
        
        lats_outer = np.where((lats < sub_lats[0]) | (lats >= sub_lats[-1]))[0]
        lons_outer = np.where((lons < sub_lons[0]) | (lons >= sub_lons[-1]))[0]
        
        all_grid[lats_outer,:] = np.nan
        all_grid[:,lons_outer] = np.nan
        
    # need to work out how to include UK DN in the current setup
    # or switch to having the mask in the format [region,lat,lon] 
    # where a new dimension is used for each mask
    #if include_uk_dn == True:
        
    #    print('Including UK devolved nations masks')
        
    # create output dataset
    if reset_index == True:
        ncountries = np.arange(len(country_names)+1)
    else:
        ncountries = [0] + indexes_all

    country_names_all = ['OCEAN'] + country_names
    country_codes_all = ['OCEAN'] + country_codes
    
    lat_da = xr.DataArray(lats,
                          coords = {'lat':lats},
                          dims = ('lat'),
                          attrs = {"long_name":"latitude","units":"degrees_north"})
    lon_da = xr.DataArray(lons,
                          coords = {'lon':lons},
                          dims = ('lon'),
                          attrs = {"long_name":"longitude","units":"degrees_east"})
    

    country_da = xr.DataArray(all_grid,
                               dims=["lat", "lon"],
                               coords=[lat_da,lon_da],
                             attrs={'long_name':'Land and ocean country indices'})

    country_names_da = xr.DataArray(country_names_all,
                                   dims=["ncountries"],
                                   coords={'ncountries':ncountries},
                                    attrs={'long_name':'Country names'})
    
    country_codes_da = xr.DataArray(country_codes_all,
                                   dims=["ncountries"],
                                   coords={'ncountries':ncountries},
                                    attrs={'long_name':'Country codes'})

    ds = xr.Dataset({"country":country_da,
                     "name":country_names_da,
                     "country_code":country_codes_da
                    })
    
    ds.attrs["Notes"] = "Created using NaturalEarth 10m Land and Marineregion datasets"
    ds.attrs["Marine_regions_data"] = "https://www.marineregions.org/eez.php"
    ds.attrs["Land_regions_data"] = "https://www.naturalearthdata.com/downloads/10m-cultural-vectors/"
    ds.attrs["Includes_ocean_territories"] = str(include_ocean_territories)
    ds.attrs["domain"] = domain
    ds.attrs["Created_by"] = f"{getpass.getuser()}@bristol.ac.uk"
    ds.attrs["Created_on"] = str(pd.Timestamp.now(tz="UTC"))
    ds.attrs['Note_on_UK_total'] = ('Isle of Man, Guernsey and Jersey have been included in the UK mask'+
                                    ' so their mask indicies will not be present in the mask')
    
    if save == True:
        if output_path is None:
            output_path = f'country_{domain}'
        
        ds.to_netcdf(f'{output_path}.nc')
        print(f'Output saved to {output_path}.')
        
    return ds           

if __name__=="__main__":
    
    ## EXAMPLE OF HOW THIS MODULE CAN BE USED ##
    write = True
    domain = "NORTHAFRICA"
    # Lat/Lon can be specified explictly or a domain footprint file can be used to find these values.
    
    ds = create_country_mask(domain,write=write,reset_index=True,ocean_label=True)
