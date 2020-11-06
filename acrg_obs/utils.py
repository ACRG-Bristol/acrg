# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:31:42 2018

@author: chxmr
"""
from __future__ import print_function

from builtins import str
from builtins import zip
import datetime
import glob
import os
import sys
from os.path import join
from datetime import datetime as dt
import json
import time
import getpass
from acrg_config.paths import paths
import xarray as xr
from pandas import Timestamp
import pandas as pd
import sqlite3

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH")
    obs_directory = os.path.join(data_path, "obs")
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg
    obs_directory = paths.obs

# Output unit strings (upper case for matching)
unit_species = {"CO2": "1e-6",
                "CH4": "1e-9",
                "C2H6": "1e-9",
                "N2O": "1e-9",
                "CO": "1e-9",
                "DCH4C13": "1",
                "DCH4D": "1",
                "DCO2C13": "1",
                "DCO2C14": "1",
                "DO2N2" : "1",
                "APO" : "1"}

# If units are non-standard, this attribute can be used
unit_species_long = {"DCH4C13": "permil",
                     "DCH4D": "permil",
                     "DCO2C13": "permil",
                     "DCO2C14": "permil",
                     "DO2N2" : "per meg",
                     "APO" : "per meg",
                     "RN": "mBq m-3"}

unit_interpret = {"ppm": "1e-6",
                  "ppb": "1e-9",
                  "ppt": "1e-12",
                  "ppq": "1e-15",
                  "else": "unknown"}

# For species which need more than just a hyphen removing and/or changing to lower case
# First element of list is the output variable name,
# second is the long name for variable standard_name and long_name
# Keys are upper case
species_translator = {"CO2": ["co2", "carbon_dioxide"],
                      "CH4": ["ch4", "methane"],
                      "ETHANE": ["c2h6", "ethane"],
                      "PROPANE": ["c3h8", "propane"],
                      "C-PROPANE": ["cc3h8", "cpropane"],
                      "BENZENE": ["c6h6", "benzene"],
                      "TOLUENE": ["c6h5ch3", "methylbenzene"],
                      "ETHENE": ["c2f4", "ethene"],
                      "ETHYNE": ["c2h2", "ethyne"],
                      "N2O": ["n2o", "nitrous_oxide"],
                      "CO": ["co", "carbon_monoxide"],
                      "H-1211": ["halon1211", "halon1211"],
                      "H-1301": ["halon1301", "halon1301"],
                      "H-2402": ["halon2402", "halon2402"],
                      "PCE": ["c2cl4", "tetrachloroethylene"],
                      "TCE": ["c2hcl3", "trichloroethylene"],
                      "PFC-116": ["c2f6", "hexafluoroethane"],
                      "PFC-218": ["c3f8", "octafluoropropane"],
                      "PFC-318": ["c4f8", "cyclooctafluorobutane"],
                      "F-113": ["cfc113", "cfc113"],
                      "H2_PDD": ["h2", "hydrogen"],
                      "NE_PDD": ["Ne", "neon"],
                      "DO2N2": ["do2n2", "ratio_of_oxygen_to_nitrogen"],
                      "DCH4C13": ["dch4c13", "delta_ch4_c13"],
                      "DCH4D": ["dch4d", "delta_ch4_d"],
                      "DCO2C13": ["dco2c13", "delta_co2_c13"],
                      "DCO2C14": ["dco2c14", "delta_co2_c14"],
                      "APO": ["apo", "atmospheric_potential_oxygen"]
                      }



def site_info_attributes(site, network):
    """ Read site information from JSON

        Args:
            site (str): Three letter site code 
            network (str): Network name
        Returns:
            dict: Dictionary containing site attributes
    """
    # Read site info file
    site_info_file = join(acrg_path, "acrg_site_info.json")
    with open(site_info_file) as sf:
        site_params = json.load(sf)

    attributes = {}
    attributes_dict = {"longitude": "station_longitude",
                       "latitude": "station_latitude",
                       "long_name": "station_long_name",
                       "height_station_masl": "station_height_masl"}

    if network is None:
        network = list(site_params[site].keys())[0]

    if site in site_params:
        for attr in attributes_dict:
            if attr in site_params[site][network]:
                attr_key = attributes_dict[attr]

                attributes[attr_key] = site_params[site][network][attr]
    else:
        # Here we return an empty dictionary as the update dictionary function will not
        # accept a NoneType
        return {}

    return attributes

def attributes(ds, species, site, network=None, global_attributes=None, units=None, scale=None, 
                                                            sampling_period=None, date_range=None):
    """ 
    This function writes attributes to an xarray.Dataset so that they conform with 
    the CF Convention v1.6
    
    Attributes of the xarray DataSet are modified, and variable names are changed

    If the species is a standard mole fraction then either:
        - species name will used in lower case in the file and variable names
            but with any hyphens taken out
        - name will be changed according to the species_translator dictionary

    If the species is isotopic data or a non-standard variable (e.g. APO):
        - Isotopes species names should begin with a "D"
            (Annoyingly, the code currently picks up "Desflurane" too. I've
             fixed this for now, but if we get a lot of other "D" species, we
             should make this better)
        - I suggest naming for isotopologues should be d<species><isotope>, e.g.
            dCH4C13, or dCO2C14
        - Any non-standard variables should be listed in the species_translator
            dictionary

    Args:
        ds (xarray.Dataset): Should contain variables such as "ch4", "ch4 repeatability".
            Must have a "time" dimension.
        species (str): Species name. e.g. "CH4", "HFC-134a", "dCH4C13"
        site (string): Three-letter site code
        network (str, default=None): Network site is associated with
        global_attribuates (dict, default=None): Dictionary containing any info you want to
            add to the file header (e.g. {"Contact": "Matt Rigby"})
        units (str, default=None): This routine will try to guess the units
            unless this is specified. Options are in units_interpret
        scale (str, default=None): Calibration scale for species. 
        sampling_period (int, default=None): Number of seconds for which air
            sample is taken. Only for time variable attribute
        date_range (list of two strings, optional): Start and end date for output
            If you only want an end date, just put a very early start date
            (e.g. ["1900-01-01", "2010-01-01"])
    """
    if not isinstance(ds, xr.Dataset):
        raise TypeError("This function only accepts xarray Datasets")

    # Current CF Conventions (v1.7) demand that valid variable names
    # begin with a letter and be composed of letters, digits and underscores
    # Here variable names are also made lowercase to enable easier matching below
    to_underscores = {var: var.lower().replace(" ", "_") for var in ds.variables}
    ds = ds.rename(to_underscores)

    species_upper = species.upper()
    species_lower = species.lower()

    matched_keys = [var for var in ds.variables if species_lower in var]

    # If we don't have any variables to rename, raise an error
    if not matched_keys:
        raise NameError(f"Cannot find species {species} in Dataset variables")
    
    species_rename = {}
    for var in matched_keys:
        if species_upper in species_translator:
            species_label = species_translator[species_upper][0]
        else:
            species_label = species_lower.replace("-", "")
        
        species_rename[var] = var.replace(species_lower, species_label)

    ds = ds.rename(species_rename)

    # Global attributes
    global_attributes_default =  {"Conditions of use": "Ensure that you contact the data owner at the outset of your project.",
                                    "Source": "In situ measurements of air",
                                    "Conventions": "CF-1.6"}

    if global_attributes:
        global_attributes.update(global_attributes_default)
    else:
        global_attributes = global_attributes_default
        
    global_attributes["File created"] = str(Timestamp.now(tz="UTC"))
    global_attributes["Processed by"] = f"{getpass.getuser()}@bristol.ac.uk"
    global_attributes["species"] = species_label

    if scale is None:
        global_attributes["Calibration_scale"] = "unknown"
    else:
        global_attributes["Calibration_scale"] = scale

    # Update the Dataset attributes
    ds.attrs.update(global_attributes)

    # Add some site attributes
    site_attributes = site_info_attributes(site.upper(), network)
    ds.attrs.update(site_attributes)

    # Species-specific attributes
    # Long name
    if (species_upper.startswith("D") and species_upper != "DESFLURANE") or species_upper == "APD":
        sp_long = species_translator[species_upper][1]
    elif species_upper == "RN":
        sp_long = "radioactivity_concentration_of_222Rn_in_air"
    elif species_upper in species_translator:
        name = species_translator[species_upper][1]
        sp_long = f"mole_fraction_of_{name}_in_air"
    else:
        sp_long = f"mole_fraction_of_{species_label}_in_air"

    ancillary_variables = []

    # Write units as attributes to variables containing any of these
    match_words = ["variability", "repeatability", "stdev", "count"]

    for key in ds.variables:
        key = key.lower()

        if species_label in key:
            # Standard name attribute
            # ds[key].attrs["standard_name"]=key.replace(species_label, sp_long)
            ds[key].attrs["long_name"] = key.replace(species_label, sp_long)

            # If units are required for variable, add attribute
            if key == species_label or any(word in key for word in match_words):
                if units is not None:
                    if units in unit_interpret:
                        ds[key].attrs["units"] = unit_interpret[units]
                    else:
                        ds[key].attrs["units"] = unit_interpret["else"]
                else:
                    ds[key].attrs["units"] = unit_species[species_upper]

                # If units are non-standard, add explanation
                if species_upper in unit_species_long:
                    ds[key].attrs["units_description"] = unit_species_long[species_upper]

            # Add to list of ancilliary variables
            if key != species_label:
                ancillary_variables.append(key)

    # TODO - for the moment skip this step - check status of ancilliary variables in standard
    # Write ancilliary variable list
    # ds[species_label].attrs["ancilliary_variables"] = ", ".join(ancillary_variables)

    # Add quality flag attributes
    # NOTE - I've removed the whitespace before status_flag and integration_flag here
    quality_flags = [key for key in ds.variables if "status_flag" in key]

    for key in quality_flags:
        ds[key] = ds[key].astype(int)
        try:
            long_name = ds[species_label].attrs["long_name"]
        except KeyError:
            raise KeyError(key, quality_flags)

        ds[key].attrs = {"flag_meaning":
                        "0 = unflagged, 1 = flagged",
                        "long_name": f"{long_name} status_flag"}

    # Add integration flag attributes
    integration_flags = [key for key in ds.variables if "integration_flag" in key]

    for key in integration_flags:
        ds[key] = ds[key].astype(int)
        long_name = ds[species_label].attrs["long_name"]
        ds[key].attrs = {"flag_meaning": "0 = area, 1 = height",
                        "standard_name": f"{long_name} integration_flag",
                        "comment": "GC peak integration method (by height or by area). Does not indicate data quality"}

    # I feel there should be a better way of doing this
    # but xarray doesn't currently have a duplicates method
    # See this https://github.com/pydata/xarray/issues/2108
    # TODO - fix this - just remove duplicates?
    if len(set(ds.time.values)) < len(ds.time.values):
        print("WARNING. Dupliate time stamps")

    first_year = Timestamp(ds.time[0].values).year

    ds.time.encoding = {"units": f"seconds since {str(first_year)}-01-01 00:00:00"}

    time_attributes = {}
    time_attributes["label"] = "left"
    time_attributes["standard_name"] = "time"
    time_attributes["comment"] = "Time stamp corresponds to beginning of sampling period. " + \
                                "Time since midnight UTC of reference date. " + \
                                "Note that sampling periods are approximate."

    if sampling_period is not None:
        time_attributes["sampling_period_seconds"] = sampling_period

    ds.time.attrs.update(time_attributes)

    # If a date range is specified, slice dataset
    if date_range:
        ds = ds.loc[dict(time = slice(*date_range))]

    return ds


def output_filename(output_directory,
                    network,
                    instrument,
                    site,
                    time_start,
                    species,
                    inlet = None,
                    version = None):
    """
    Create an output filename in the format
    output_directory/site/network-instrument_site_YYYYMMDD_species[-inlet]-version.nc

    Also creates a site directory in output_directory, if one doesn't exist

    Args:

        inlet (string, optional): Label for inlet. If none supplied, assumes that there
            is only one inlet, and no inlet info is added to file name
        version (string, optional): Version number for file. If none supplied, will
            use YYYYMMDD for day processed

    """

    # Check if directory exists. If not, create it
    if not os.path.exists(join(output_directory, site)):
        os.makedirs(join(output_directory, site))

    # Create suffix. If inlet is specified, append to species name
    suffix = species
    if inlet is not None:
        suffix += "-" + inlet
    if version is not None:
        suffix += "-" + version
    else:
        suffix += "-" + time.strftime("%Y%m%d")

    # Return output filename
    return join(output_directory,
                "%s/%s-%s_%s_%04i%02i%02i_%s.nc" %(site,
                                                   network,
                                                   instrument,
                                                   site,
                                                   time_start.year,
                                                   time_start.month,
                                                   time_start.day,
                                                   suffix))


def obs_database(data_directory = None):
    '''
    Creates an SQLite database in obs folder detailing contents of each file
    
    Args:
        data_directory (pathlib, optional) : Path to user-defined obs_folder
    '''
    
    # Directories to exclude from database
    exclude = ["GOSAT", "unknown"]

    network = []
    instrument = []
    site_code = []
    start_date = []
    end_date = []
    species = []
    inlet = []
    calibration_scale = []
    filename = []

    if data_directory is None:
        obs_path = paths.obs
    else:
        obs_path = data_directory
    
    print(f"Reading obs files in {obs_path}")
    
    # Find sub-directories in obs folder
    for d in obs_path.glob("*"):
        if d.is_dir():
            if d.name not in exclude:

                # Find netcdf files
                files = d.glob("*.nc")
                for f in files:

                    # TODO: Some files are empty. Figure out why!
                    if os.stat(f).st_size != 0:

                        #TODO: add try/except to see if files open with xarray 

                        f_parts = f.name.split("_")

                        network.append(f_parts[0].split("-")[0])
                        instrument.append(f_parts[0].split("-")[1])

                        site_code.append(f_parts[1])

                        extras = f_parts[3].split("-")
                        species.append(extras[0])
                        if len(extras) == 3:
                            inlet.append(extras[1])
                        else:
                            inlet.append("%")

                        with xr.open_dataset(f) as ds:
                            if "Calibration_scale" in ds.attrs.keys():
                                calibration_scale.append(ds.attrs["Calibration_scale"])
                            else:
                                calibration_scale.append(None)
                            start_date.append(pd.Timestamp(ds["time"].values[0]).to_pydatetime())
                            end_date.append(pd.Timestamp(ds["time"].values[-1]).to_pydatetime())

                        filename.append(str(f))
                        
    file_info = list(zip(filename, network, instrument, site_code, species, inlet, calibration_scale, start_date, end_date))    

    # Write database
    ######################################
    
    print(f"Writing database {obs_path / 'obs.db'}")
    
    conn = sqlite3.connect(obs_path / "obs.db")
    c = conn.cursor()

    c.execute('''DROP TABLE IF EXISTS files
              ''')
    
    c.execute('''CREATE TABLE files
                 (filename text, network text, instrument text, site text, species text, inlet text, scale text, startDate timestamp, endDate timestamp)''')

    c.executemany('INSERT INTO files VALUES (?,?,?,?,?,?,?,?,?)', file_info)

    conn.commit()
    conn.close()



def cleanup(site,
            version = None):
    '''
    Archive old versions of files in an archive folder within each site folder.
    Will keep the maxiumum version number for each NETWORK-INSTRUMENT, and move
    older ones to a zip file.

    Args:
        site (string): Measurement site

        version (string): Specify a version to archive. If None, will archive older versions
    '''


    import zipfile

    # Directories and files
    data_directory = os.path.join(obs_directory, site)
    archive_directory = os.path.join(data_directory, "archive")
    files = glob.glob(os.path.join(data_directory, "*.nc"))

    print("Checking for files to archive in %s ..." % data_directory)

    # Work out instrument (NETWORK-INSTRUMENT) and suffix from filename
    instrument = []
    suffix = []

    for f in files:
        f_data = os.path.split(f)[1].split(".")[0].split("_")
        instrument.append(f_data[0])
        suffix.append(f_data[3])

    # Identify unique instruments
    unique_instruments = set(instrument)

    # Cycle through each NETWORK-INSTRUMENT
    for i in unique_instruments:

        # Find versions
        versions = [s.split("-")[-1] for (s, inst) in zip(suffix, instrument) if inst == i]
        unique_versions = set(versions)
        latest_version = max(unique_versions)

        # By default, archive everything apart from the latest version
        archive_versions = [v for v in unique_versions if v != latest_version]

        # If a particular version is specified, only archive that one
        if version is not None:
            if version not in unique_versions:
                print("Error: version %s not in %s. Versions in this folder: %s" %
                      (version, data_directory, ",".join(unique_versions)))
                return
            archive_versions = [version]

        if len(archive_versions) == 0:
            print("... everything up-to-date for %s" % i)

        # For any old versions, archive and delete files
        for v in archive_versions:

            print("... archiving %s files, version %s" % (i, v) )

            # Check if archive directory exists, and create if not
            if not os.path.exists(archive_directory):
                os.makedirs(archive_directory)

            # Archive filepath for this instrument and version
            archive = os.path.join(archive_directory,
                                   "%s_%s_%s.zip" % (i, site, v))

            # Files to archive for this instrument and version
            archive_files = [f for f in files if "-" + v in f and i + "_" in f]

            # Write archive and delete file
            zipf = zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED)
            for f in archive_files:
                zipf.write(f, os.path.basename(f))
                os.remove(f)
            zipf.close()

            
def cleanup_all():
    
    site_list = [site_dir for site_dir in os.listdir(obs_directory) if os.path.isdir(os.path.join(obs_directory,site_dir))]
    for site in site_list:
        cleanup(site)
