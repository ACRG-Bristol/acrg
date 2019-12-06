#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:31:42 2018

@author: chxmr
"""
from __future__ import print_function

from builtins import str
from builtins import zip
import glob
import os
from os.path import join
from datetime import datetime as dt
import json
import time
import getpass

#acrg_path = os.path.dirname(os.path.realpath(__file__))
acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

# Set default obs folder
obs_directory = os.path.join(data_path, "obs/")

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

    # Read site info file
    site_info_file = join(acrg_path, "acrg_site_info.json")
    with open(site_info_file) as sf:
        site_params = json.load(sf)

    attributes = {}
    attributes_list = {"longitude": "station_longitude",
                       "latitude": "station_latitude",
                       "long_name": "station_long_name",
                       "height_station_masl": "station_height_masl"}

    if network is None:
        net = list(site_params[site].keys())[0]
    else:
        net = network

    if site in list(site_params.keys()):
        for at in list(attributes_list.keys()):
            if at in list(site_params[site][net].keys()):
                attributes[attributes_list[at]] = site_params[site][net][at]
        return attributes
    else:
        return None

def attributes(ds, species, site,
               network = None,
               global_attributes = None,
               units = None,
               scale = None,
               sampling_period = None,
               date_range = None,
               global_attributes_default = {"Conditions of use": "Ensure that you contact the data owner at the outset of your project.",
                                            "Source": "In situ measurements of air",
                                            "Conventions": "CF-1.6"}):
    """
    Format attributes for netCDF file
    Attributes of xarray DataSet are modified, and variable names are changed

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
        ds (xarray dataset): Should contain variables such as "ch4", "ch4 repeatability".
            Must have a "time" dimension.
        species (string): Species name. e.g. "CH4", "HFC-134a", "dCH4C13"
        site (string): Three-letter site code

        global_attribuates (dict, optional): Dictionary containing any info you want to
            add to the file header (e.g. {"Contact": "Matt Rigby"})
        units (string, optional): This routine will try to guess the units
            unless this is specified. Options are in units_interpret
        scale (string, optional): Calibration scale for file header.
        sampling_period (int, optional): Number of seconds for which air
            sample is taken. Only for time variable attribute
        date_range (list of two strings, optional): Start and end date for output
            If you only want an end date, just put a very early start date
            (e.g. ["1900-01-01", "2010-01-01"])
    """

    # Rename all columns to lower case! Could this cause problems?
    rename_dict = {}
    for key in ds.variables:
        rename_dict[key] = key.lower()
    ds = ds.rename(rename_dict)

    # Rename species, if required
    rename_dict = {}
    for key in ds.variables:
        if species.lower() in key:
            if species.upper() in list(species_translator.keys()):
                # Rename based on species_translator, if available
                species_out = species_translator[species.upper()][0]
            else:
                # Rename species to be lower case and without hyphens
                species_out = species.lower().replace("-", "")

            rename_dict[key] = key.replace(species.lower(), species_out)
    ds = ds.rename(rename_dict)

    # Check if these was a variable with the species name in it
    try:
      species_out
    except NameError:
      print("ERROR: Can't find species %s in column names %s" %(species, list(ds.keys())))

    # Global attributes
    #############################################

    if global_attributes is None:
        global_attributes = {}
        for key in global_attributes_default:
            global_attributes[key] = global_attributes_default[key]
    else:
        for key in global_attributes_default:
            global_attributes[key] = global_attributes_default[key]

    # Add some defaults
    global_attributes["File created"] = str(dt.now())

    # Add user
    global_attributes["Processed by"] = "%s@bristol.ac.uk" % getpass.getuser()


    for key, values in global_attributes.items():
        ds.attrs[key] = values

    # Add some site attributes
    global_attributes_site = site_info_attributes(site.upper(), network)
    if global_attributes_site is not None:
        for key, values in global_attributes_site.items():
            ds.attrs[key] = values

    # Add calibration scale
    if scale:
        ds.attrs["Calibration_scale"] = scale
    else:
        ds.attrs["Calibration_scale"] = "unknown"

    # Add species name
    ds.attrs["species"] = species_out

    # Species-specific attributes
    #############################################

    # Long name
    if (species.upper()[0] == "D" and species.upper() != "DESFLURANE") or species.upper() == "APO":
        sp_long = species_translator[species.upper()][1]
    elif species.upper() == "RN":
        sp_long = "radioactivity_concentration_of_222Rn_in_air"
    elif species.upper() in list(species_translator.keys()):
        sp_long = "mole_fraction_of_" + species_translator[species.upper()][1] + "_in_air"
    else:
        sp_long = "mole_fraction_of_" + species_out + "_in_air"

    ancillary_variables = ""

    for key in ds.variables:

        if species_out in key:

            # Standard name attribute
            #ds[key].attrs["standard_name"]=key.replace(species_out, sp_long)
            ds[key].attrs["long_name"]=key.replace(species_out, sp_long)

            # If units are required for variable, add attribute
            if (key == species_out) or \
                ("variability" in key) or \
                ("repeatability" in key):
                if units is None:
                    ds[key].attrs["units"] = unit_species[species.upper()]
                else:
                    if units in list(unit_interpret.keys()):
                        ds[key].attrs["units"] = unit_interpret[units]
                    else:
                        ds[key].attrs["units"] = unit_interpret["else"]

                # if units are non-standard, add explanation
                if species.upper() in list(unit_species_long.keys()):
                    ds[key].attrs["units_description"] = unit_species_long[species.upper()]

            # Add to list of ancilliary variables
            if key != species_out:
                ancillary_variables += key + ", "

    # Write ancilliary variable list
    ds[species_out].attrs["ancilliary_variables"] = ancillary_variables.strip()

    # Add quality flag attributes
    ##################################

    flag_key = [key for key in ds.variables if " status_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = unflagged, 1 = flagged",
                              "long_name":
                              ds[species_out].attrs["long_name"] + " status_flag"}

    # Add integration flag attributes
    ##################################

    flag_key = [key for key in ds.variables if " integration_flag" in key]
    if len(flag_key) > 0:
        flag_key = flag_key[0]
        ds[flag_key] = ds[flag_key].astype(int)
        ds[flag_key].attrs = {"flag_meaning":
                              "0 = area, 1 = height",
                              "standard_name":
                              ds[species_out].attrs["long_name"] + " integration_flag",
                              "comment":
                              "GC peak integration method (by height or by area). " +
                              "Does not indicate data quality"}

    # Set time encoding
    #########################################

    # Check if there are duplicate time stamps
    if len(set(ds.time.values)) < len(ds.time.values):
        print("WARNING. Dupliate time stamps")

    first_year = str(ds.time.to_pandas().index.to_pydatetime()[0].year)

    ds.time.encoding = {"units": "seconds since " + \
                        first_year + "-01-01 00:00:00"}
    ds.time.attrs["label"] = "left"
    ds.time.attrs["comment"] = "Time stamp corresponds to beginning of sampling period. " + \
                               "Time since midnight UTC of reference date. " + \
                               "Note that sampling periods are approximate."
    if sampling_period:
        ds.time.attrs["sampling_period_seconds"] = sampling_period

    # If a date range is specified, slice dataset
    if date_range != None:
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
