# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 14:08:07 2015

@author: chxmr
"""
import datetime
import numpy as np
import pandas as pd
from os.path import join, split
from datetime import datetime as dt
import glob
import xarray as xray
import json
from os import stat
import fnmatch
from .utils import attributes, output_filename, cleanup

from acrg.config.paths import Paths

acrg_path = Paths.acrg

# Read site info file
info_file = join(acrg_path,
                 "acrg/obs/process_gcwerks_parameters.json")
with open(info_file) as sf:
    params = json.load(sf)

site_info_file = join(acrg_path, "data/site_info.json")
with open(site_info_file) as sf:
    site_params = json.load(sf)

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
                     "APO" : "per meg"}

unit_interpret = {"ppm": "1e-6",
                  "ppb": "1e-9",
                  "ppt": "1e-12",
                  "ppq": "1e-15",
                  "else": "unknown"}


# For species which need more than just a hyphen removing or changing to lower case
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

def parser_YYMMDD(yymmdd):
    return dt.strptime(yymmdd, '%y%m%d')


def get_directories(default_input_directory,
                    default_output_directory,
                    user_specified_input_directory = None,
                    user_specified_output_directory = None):

    # If an output directory directory is set, use that, otherwise from json file
    if user_specified_output_directory:
        output_folder = user_specified_output_directory
    else:
        output_folder = default_output_directory

    # If an input directory directory is set, use that, otherwise from json file
    if user_specified_input_directory:
        if user_specified_output_directory is None:
            print("EXITING: You must also specify an output folder if you specify an input")
            return None
        data_folder = user_specified_input_directory
    else:
        data_folder = default_input_directory

    return data_folder, output_folder


# ICOS
########################################################

def icos_data_read(data_file, species):

    print("Reading " + data_file)

    # Find out how many header lines there are
    nheader = 0
    with open(data_file, "r") as f:
        for l in f:
            if l[0] != "#":
                break
            nheader += 1

    # Read CSV file
    df =  pd.read_csv(data_file,
                      skiprows = nheader-1,
                      parse_dates = {"time": ["Year", "Month", "Day", "Hour", "Minute"]},
                      index_col = "time",
                      sep = ";",
                      usecols = ["Day", "Month", "Year", "Hour", "Minute",
                                 str(species.lower()), "SamplingHeight",
                                 "Stdev", "NbPoints"],
                      dtype = {"Day": np.int,
                               "Month": np.int,
                               "Year": np.int,
                               "Hour": np.int,
                               "Minute": np.int,
                               species.lower(): np.float,
                               "Stdev": np.float,
                               "SamplingHeight": np.float,
                               "NbPoints": np.int},
                      na_values = "-999.99")

    # Format time
    df.index = pd.to_datetime(df.index, format = "%Y %m %d %H %M")

    df = df[df[species.lower()] >= 0.]

    # Remove duplicate indices
    df.reset_index(inplace = True)
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')

    # Rename columns
    df.rename(columns = {species.lower(): species.upper(),
                         "Stdev": species.upper() + " variability",
                         "NbPoints": species.upper() + " number_of_observations"},
               inplace = True)

    df.index.name = "time"

    # Convert to Dataset
    ds = xray.Dataset.from_dataframe(df.sort_index())

    return ds


def icos(site, network = "ICOS",
         input_directory = None,
         output_directory = None,
         date_range = None,
         version = None):

    def find_species_inlet_model(filenames):
        out = []
        for f in filenames:
            f_elements = f.split(".")
            if len(f_elements) == 6:
                out.append((f_elements[1],
                            f_elements[4],
                            "picarro" + f_elements[3].upper()))
            else:
                out.append((f_elements[1],
                            f_elements[3],
                            "picarro"))
        return(out)

    # Get directories and site strings
    params_icos = params["ICOS"]
    site_string = params_icos[site]["gcwerks_site_name"]

    data_folder, output_folder = \
            get_directories(params_icos["directory"].replace("%site", site_string),
                            params_icos["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    # Search for species, inlets and model from file names
    data_file_search = join(data_folder, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    data_file_names = [split(f)[1] for f in data_files]
    species_inlet_model = find_species_inlet_model(data_file_names)

    inlets = set([i for (s, i, m) in species_inlet_model])

    for i, (species, inlet, model) in enumerate(species_inlet_model):

        if stat(data_files[i]).st_size > 0:

            # Create Pandas dataframe
            ds = icos_data_read(data_files[i], species.upper())

            # Sort out attributes
            global_attributes = params_icos[site.upper()]["global_attributes"]
            global_attributes["inlet_height_magl"] = float(params_icos[site]["inlet_rename"][inlet][:-1])

            ds = attributes(ds,
                            species.upper(),
                            site.upper(),
                            network = network,
                            global_attributes = global_attributes,
                            sampling_period = 60,
                            date_range = date_range)

            if len(ds.time.values) == 0:

                # Then must have not passed date_range filter?
                print(" ... no data in range")
                # then do nothing

            else:

                inlet_label = params_icos[site]["inlet_rename"][inlet]

                # Write file
                nc_filename = output_filename(output_folder,
                                              network,
                                              model,
                                              site.upper(),
                                              ds.time.to_pandas().index.to_pydatetime()[0],
                                              ds.species,
                                              inlet = [None, inlet_label][len(inlets) > 1],
                                              version = version)

                ds.to_netcdf(nc_filename)

                print("Written " + nc_filename)

        else:
            print("Skipping empty file: %s" % data_files[i])

# GC FUNCTIONS
###############################################################

def gc_data_read(dotC_file, scale = {}, units = {}):
    """ Read GC data

        Args:
            dotC_file (str): Filepath for .C data file
            scale (dict): Dictionary to hold scales for each species
            units (dict): Dictionary to hold units for each species
        Returns:
            tuple: Tuple of pandas.Dataframe, species list, units dict, species dict
    """
    # Read header
    header = pd.read_csv(dotC_file,
                         skiprows=2,
                         nrows=2,
                         header = None,
                         sep=r"\s+")

     # Create a function to parse the datetime in the data file
    def parser(date): 
        return datetime.datetime.strptime(date, '%Y %m %d %H %M')
    # Read the data in and automatically create a datetime column from the 5 columns
    # Not dropping the yyyy', 'mm', 'dd', 'hh', 'mi' columns here
    df = pd.read_csv(dotC_file, 
                    skiprows=4, 
                    sep=r"\s+", 
                    index_col=["yyyy_mm_dd_hh_mi"],
                    parse_dates=[[1, 2, 3, 4, 5]],
                    date_parser=parser, 
                    keep_date_col=True)

    df.index.name = "time"
    # Drop duplicates
    df = df.reset_index().drop_duplicates(subset='time').set_index('time')

    units = {}
    scale = {}
    species = []

    columns_renamed = {}
    for column in df.columns:
        if "Flag" in column:
            # Location of this column in a range (0, n_columns-1)
            col_loc = df.columns.get_loc(column)
            # Get name of column before this one for the gas name
            gas_name = df.columns[col_loc - 1]
            # Add it to the dictionary for renaming later
            columns_renamed[column] = gas_name + "_flag"
            # Create 2 new columns based on the flag columns
            df[gas_name + " status_flag"] = (df[column].str[0] != "-").astype(int)
            df[gas_name + " integration_flag"] = (df[column].str[1] != "-").astype(int)

            col_shift = -1
            units[gas_name] = header.iloc[1, col_loc + col_shift]
            scale[gas_name] = header.iloc[0, col_loc + col_shift]

            # Ensure the units and scale have been read in correctly
            # Have this in case the column shift between the header and data changes
            if units[gas_name] == "--" or scale[gas_name] == "--":
                raise ValueError("Error reading units and scale, ensure columns are correct between header and dataframe")

            species.append(gas_name)

    # Rename columns to include the gas this flag represents
    df = df.rename(columns=columns_renamed, inplace=False)

    return df, species, units, scale

def gc_precisions_read(precisions_file):

    # Read precision species
    precision_species = list(pd.read_csv(precisions_file,
                                         skiprows=3,
                                         nrows = 1,
                                         header = None,
                                         sep=r"\s+").values[0][1:])

    # Read precisions
    precision = pd.read_csv(precisions_file,
                            skiprows=5,
                            header = None,
                            sep=r"\s+", dtype = str,
                            index_col = 0,
                            parse_dates = True,
                            date_parser = parser_YYMMDD)

    # Rename index column
    precision.index.names = ["index"]

    # Drop duplicates
    precision = precision.reset_index().drop_duplicates(subset='index').set_index('index')

    return precision, precision_species


def gc(site, instrument, network,
       input_directory = None,
       output_directory = None,
       date_range = None,
       version = None):
    """
    Process GC data per site and instrument
    Instruments can be:
        "GCMD": GC multi-detector (output will be labeled GC-FID or GC-ECD)
        "GCMS": GC ADS (output GC-ADS)
        "medusa": GC medusa (output GC-MEDUSA)

    Network is the network name for output file.
    """

    site_gcwerks = params["GC"][site]["gcwerks_site_name"]
    instrument_gcwerks = params["GC"]["instruments"][instrument]

    data_folder, output_folder = \
            get_directories(params["GC"]["directory"][instrument],
                            params["GC"]["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    search_strings = []
    for suffix in params["GC"]["instruments_suffix"][instrument]:
        # Search string
        search_string = join(data_folder,
                             site_gcwerks + \
                             instrument_gcwerks + \
                             suffix + ".??.C")
        search_strings.append(search_string)

        data_files = sorted(glob.glob(search_string))
        if len(data_files) > 0:
            break

    # Error if can't find files
    if len(data_files) == 0.:
        print("ERROR: can't find any files: " + \
              ",\r".join(search_strings))
        return None

    precision_files = [data_file[0:-2] + ".precisions.C" \
                        for data_file in data_files]

    dfs = []
    scale = {}
    units = {}

    for fi, data_file in enumerate(data_files):

        print("Reading " + data_file)

        # Get observations
        df, species, units, scale = gc_data_read(data_file,
                                                 scale = scale,
                                                 units = units)

        # Get precision
        precision, precision_species = gc_precisions_read(precision_files[fi])

        # Merge precisions into dataframe
        for sp in species:
            precision_index = precision_species.index(sp)*2+1
            df[sp + " repeatability"] = precision[precision_index].\
                                            astype(float).\
                                            reindex_like(df, "pad")
        dfs.append(df)

    # Concatenate
    dfs = pd.concat(dfs).sort_index()

    # Apply timestamp correction, because GCwerks currently outputs
    #   the CENTRE of the sampling period
    dfs["new_time"] = dfs.index - \
            pd.Timedelta(seconds = params["GC"]["sampling_period"][instrument]/2.)
    dfs = dfs.set_index("new_time", inplace=False, drop=True)

    # Label time index
    dfs.index.name = "time"

    # Convert to xray dataset
    ds = xray.Dataset.from_dataframe(dfs)

    # Get species from scale dictionary
    species = list(scale.keys())

    inlets = params["GC"][site]["inlets"]

    # Process each species in file
    for sp in species:

        global_attributes = params["GC"][site.upper()]["global_attributes"]
        global_attributes["comment"] = params["GC"]["comment"][instrument]

        # Now go through each inlet (if required)
        for inleti, inlet in enumerate(inlets):

            # There is only one inlet, just use all data, and don't lable inlet in filename
            if (inlet == "any") or (inlet == "air"):

                print("Processing %s, assuming single inlet..." %sp)

                ds_sp = ds[[sp,
                            sp + " repeatability",
                            sp + " status_flag",
                            sp + " integration_flag",
                            "Inlet"]]

                # No inlet label in file name
                inlet_label = None

            else:
                # Get specific inlet

                print("Processing " + sp + ", " + inlet + "...")

                # if inlet is in the format "date_YYYYMMDD_YYYYMMDD", split by date
                if inlet[0:4] == "date":
                    dates = inlet.split("_")[1:]
                    slice_dict = dict(time = slice(dates[0], dates[1]))
                    ds_sliced = ds.loc[slice_dict]
                    ds_sp = ds_sliced[[sp,
                                       sp + " repeatability",
                                       sp + " status_flag",
                                       sp + " integration_flag",
                                       "Inlet"]]

                else:

                    # Use UNIX pattern matching to find matching inlets
                    # select_inlet is a list of True or False
                    select_inlet = [fnmatch.fnmatch(i, inlet) for i in ds.Inlet.values]
                    # now create a DataArray of True or False
                    select_ds = xray.DataArray(select_inlet, coords = [ds.time],
                                               dims = ["time"])

                    # sub-set ds
                    ds_sp = ds.where(select_ds, drop = True)[[sp,
                                                              sp + " repeatability",
                                                              sp + " status_flag",
                                                              sp + " integration_flag",
                                                              "Inlet"]]

                # re-label inlet if required
                if "inlet_label" in list(params["GC"][site].keys()):
                    inlet_label = params["GC"][site]["inlet_label"][inleti]
                else:
                    inlet_label = inlet

            if inlet_label == None:
                global_attributes["inlet_magl"] = params["GC"][site]["inlet_label"][inleti]
            else:
                global_attributes["inlet_magl"] = inlet_label

            # Record Inlets from the .C file, for the record
            # TODO: figure out why xarray raises an error at this line
            #   if "analysis time" column is included (commented out above)
            Inlets = set(ds_sp.where(ds_sp[sp + " status_flag"] == 0, drop = True).Inlet.values)
            global_attributes["inlet_gcwerks"] = ", ".join(Inlets)
            # Now remove "Inlet" column from dataframe. Don't need it
            ds_sp = ds_sp.drop(["Inlet"])

            # Drop NaNs
            ds_sp = ds_sp.dropna("time")

            if len(ds_sp.time) == 0:

                print("... no data in file, skipping " + sp)

            else:

                # Sort out attributes
                ds_sp = attributes(ds_sp, sp, site.upper(),
                                   network = network,
                                   global_attributes = global_attributes,
                                   units = units[sp],
                                   scale = scale[sp],
                                   sampling_period = params["GC"]["sampling_period"][instrument],
                                   date_range = date_range)

                if len(ds_sp.time.values) == 0:

                    # Then must have not passed date_range filter?
                    print(" ... no data in range")
                    # then do nothing

                else:

                    # Get instrument name for output
                    if sp.upper() in params["GC"]["instruments_out"][instrument]:
                        instrument_out = params["GC"]["instruments_out"][instrument][sp]
                    else:
                        instrument_out = params["GC"]["instruments_out"][instrument]["else"]

                    # Write file
                    nc_filename = output_filename(output_folder,
                                                  network,
                                                  instrument_out,
                                                  site.upper(),
                                                  ds_sp.time.to_pandas().index.to_pydatetime()[0],
                                                  ds_sp.species,
                                                  inlet = [None, inlet_label][len(inlets) > 1],
                                                  version = version)

                    # compress
                    ds_sp = set_encoding(ds_sp)
                    
                    print("Writing... " + nc_filename)
                    ds_sp.to_netcdf(nc_filename)
                    print("... written.")



def crds_data_read(data_file):
    """ Read CRDS data

        Args:
            data_file (str): Path to data file
        Returns:
            tuple: Tuple of xarray.Dataset and list of species
    """

    # Translate header strings
    crds_header_string_interpret = {"C": "",
                                    "stdev": "_variability",
                                    "N": "_number_of_observations"}
    
    print("Reading " + data_file)

    # Read file header
    df_header = pd.read_csv(data_file,
                         skiprows=1,
                         nrows = 2,
                         header = None,
                         sep=r"\s+")

    header = []
    species = []

    # Create header list
    for i in df_header.columns:
        if df_header[i][0] != '-':
            header.append(df_header[i][0].upper() + \
                          crds_header_string_interpret[df_header[i][1]])
            if df_header[i][1] == "C":
                species.append(df_header[i][0].upper())
        else:
            header.append(df_header[i][1].upper())

    # Function to parse the datetime format found in the datafile
    def parse_date(date):
        try:
            return datetime.datetime.strptime(date, '%y%m%d %H%M%S')
        except ValueError:
            return pd.NaT

    # Read data
    df = pd.read_csv(data_file,
                     skiprows=4,
                     header = None,
                     sep=r"\s+",
                     names = header,
                     index_col=["DATE_TIME"],
                     parse_dates=[["DATE", "TIME"]],
                     date_parser=parse_date)

    # Check if the index is sorted and if not sort it
    if not df.index.is_monotonic_increasing:
        df.sort_index()

    df.index.name = "time"

    # Remove duplicates
    df = df.reset_index().drop_duplicates(subset='time').set_index('time')

    # Convert to a Dataset
    ds = df.to_xarray()

    return ds, species


def crds(site, network,
         input_directory = None,
         output_directory = None,
         date_range = None,
         version = None):
    """
    Process CRDS data

    Args:
        site (str): Three letter site code
        network (str): Network string only for output
        input_directory (str, default=None): Directory holding input files
        output_directory (str, default=None): Directory for output files
        date_range (list, default=None): Start and end date for output. If you only want an end date, 
                                        just put a very early start date (e.g. ["1900-01-01", "2010-01-01"])
        version (str, default=None): Version string for output filename
    """
    params_crds = params["CRDS"]

    site_string = params_crds[site]["gcwerks_site_name"]

    data_folder, output_folder = \
            get_directories(params_crds["directory"].replace("%site", site_string),
                            params["CRDS"]["directory_output"],
                            user_specified_input_directory = input_directory,
                            user_specified_output_directory = output_directory)

    # Search for species and inlets from file names
    data_file_search = join(data_folder, site.lower() + ".*.1minute.*.dat")
    data_files = glob.glob(data_file_search)
    instruments = [f.split(".")[1] for f in data_files]
    inlets = [f.split(".")[-2] for f in data_files]

    for i, inlet in enumerate(inlets):

        # Create Pandas dataframe
        ds, species = crds_data_read(data_files[i])

        # Write netCDF file for each species
        for sp in species:

            # Species-specific dataset
            ds_sp = ds[[sp,
                        sp + "_variability",
                        sp + "_number_of_observations"]]
            ds_sp = ds_sp.dropna("time")

            global_attributes = params_crds[site]["global_attributes"]
            global_attributes["inlet_height_magl"] = float(inlet[0:-1])
            global_attributes["comment"] = params_crds["comment"]

            ds_sp = attributes(ds_sp, sp, site.upper(),
                               network = network,
                               global_attributes = global_attributes,
                               scale = params_crds["default_scales"][sp],
                               sampling_period=60,
                               date_range = date_range)

            if len(ds_sp.time.values) == 0:

                # Then must have not passed date_range filter?
                print(" ... no data in range")
                # then do nothing

            else:

                # Write file
                nc_filename = output_filename(output_folder,
                                              network,
                                              instruments[i],
                                              site.upper(),
                                              ds_sp.time.to_pandas().index.to_pydatetime()[0],
                                              ds_sp.species,
                                              inlet = [None, inlet][len(inlets) > 1],
                                              version = version)
                
                # compress data
                ds_sp = set_encoding(ds_sp)
                
                print("Writing " + nc_filename)
                ds_sp.to_netcdf(nc_filename)
                print("... written.")

def set_encoding(ds):
    '''
    Specify encoding to prevent files getting too big
    
    This function makes sure that number_of_observations is an integer
    non-time variables are float32, and applies compression.
    '''
    
    for var in ds:
        if "number_of_observations" in var:
            ds[var].values = ds[var].values.astype(int)
            ds[var].encoding["dtype"]="int16"
        else:
            if var != "time":
                ds[var].encoding["dtype"] = "float32"
#            ds[var].encoding["zlib"] = True
        ds[var].encoding["zlib"]=True

    return(ds)

def ale_gage(site, network):
    """
    Process Georgia Tech ALE or GAGE observations

    Args:
        site (str): ADR, RPB, ORG, SMO, CGO or MHD
        network (str): ALE or GAGE

    """

    import fortranformat as ff

    ale_directory = "/dagage2/agage/summary/git/ale_new/complete/"
    gage_directory = "/dagage2/agage/summary/git/gage_new/complete/"

    output_directory = "/data/shared/obs_2018/"

    site_translate = {"ADR": "adrigole",
                      "RPB": "barbados",
                      "ORG": "oregon",
                      "SMO": "samoa",
                      "CGO": "tasmania",
                      "MHD": "macehead"}

    if network == "ALE":
        data_directory = ale_directory
    elif network == "GAGE":
        data_directory = gage_directory
    else:
        print("Network needs to be ALE or GAGE")
        return None


    fnames = sorted(glob.glob(join(data_directory,
                                   site_translate[site] + "/" + site + "*.dap")))


    formatter = ff.FortranRecordReader('(F10.5, 2I4,I6, 2I4,I6,1X,10(F10.3,a1))')


    dfs = []
    for fname in fnames:

        print("Reading... " + fname)

        header = []

        with open(fname) as f:
            for i in range(6):
                header.append(f.readline())

            lines = f.readlines()

        scales = header[-3].split()
        units = header[-2].split()
        species = header[-1].split()

        dayi = species.index("DD")
        monthi = species.index("MM")
        yeari = species.index("YYYY")
        houri = species.index("hh")
        mini = species.index("min")

        data = []
        time = []

        for line in lines:
            data_line = formatter.read(line)

            if data_line[mini] < 60 and data_line[houri] < 24:
                data.append([d for d in data_line if d != " " and \
                                                     d != None and \
                                                     d != "P"])
                time.append(dt(data_line[yeari],
                               data_line[monthi],
                               data_line[dayi],
                               data_line[houri],
                               data_line[mini]))

        data = np.vstack(data)
        data = data[:, 7:]

        df = pd.DataFrame(data = data, columns = species[7:], index = time)
        df.replace(to_replace = 0., value=np.NaN, inplace = True)
        dfs.append(df)

    scales = scales[7:]
    units = units[7:]
    species = species[7:]

    df = pd.concat(dfs)

    # Write netCDF file for each species
    for si, sp in enumerate(species):

        # Remove duplicate indices
        df_sp = df[sp].reset_index().drop_duplicates(subset='index').set_index('index')

        # Convert to Dataset
        df_sp.index.name = "time"
        ds = xray.Dataset.from_dataframe(df_sp.sort_index())
        ds = ds.dropna("time")

        ds = attributes(ds, sp, site.upper(),
                        network = network,
                       scale = scales[si],
                       sampling_period=60,
                       units = units[si])

        # Write file
        nc_filename = output_filename(output_directory,
                                      network,
                                      "GCECD",
                                      site.upper(),
                                      ds.time.to_pandas().index.to_pydatetime()[0],
                                      ds.species)
        print("Writing " + nc_filename)
        ds.to_netcdf(nc_filename)
        print("... written.")


def mhd_o3():

    channels = ["channel1", "channel0", "channel2"]
    base_directory = "/dagage2/agage/macehead-ozone/results/export/"

    df = []

    for channel in channels:

        files_channel = sorted(glob.glob(join(base_directory, channel, "*.csv")))

        for f in files_channel:

            df.append(pd.read_csv(f, sep=",",
                                  names = ["datetime",
                                           "ozone",
                                           "ozone_variability",
                                           "ozone_number_samples"],
                                  na_values = "NA",
                                  index_col = "datetime",
                                  parse_dates = ["datetime"]))
            df[-1].dropna(inplace = True)

    df = pd.concat(df)
    df.index.name = "index"
    df = df.reset_index().drop_duplicates(subset='index').set_index('index')
    df.sort_index(inplace = True)

    # Convert to Dataset
    df.index.name = "time"
    ds = xray.Dataset.from_dataframe(df)

    ds = attributes(ds,
                    "ozone",
                    "MHD",
                    network = network,
                    scale = "SCALE",
                    sampling_period=60*60,
                    units = "ppb")

    # Write file
    nc_filename = output_filename("/dagage2/agage/metoffice/processed_observations_2018",
                                  "AURN",
                                  "thermo",
                                  "MHD",
                                  str(ds.time.to_pandas().index.to_pydatetime()[0].year),
                                  ds.species,
                                  site_params["MHD"]["AGAGE"]["height"][0])
    print("Writing " + nc_filename)
    ds.to_netcdf(nc_filename)
    print("... written.")


def data_freeze(version,
                end_date,
                input_directory = "/dagage2/agage/summary/gccompare-net/snapshot/current-frozendata/data-net/",
                output_directory = "/dagage2/agage/summary/gccompare-net/snapshot/current-frozendata/data-net/processed/"):

    date_range = ["19000101", end_date]

    # ICOS
    icos("MHD", network = "ICOS", input_directory = input_directory,
         output_directory = output_directory, version = version, date_range = date_range)

    # GAUGE CRDS data
    crds("HFD", "GAUGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    crds("BSD", "GAUGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)


    # GAUGE GC data
    gc("BSD", "GCMD", "GAUGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)
    gc("HFD", "GCMD", "GAUGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # DECC CRDS data
    crds("TTA", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)
    crds("RGL", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)
    crds("TAC", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # DECC GC data
    gc("TAC", "GCMD", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)
    gc("RGL", "GCMD", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # DECC Medusa
    gc("TAC", "medusa", "DECC", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # AGAGE GC data
    gc("MHD", "GCMD", "AGAGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # AGAGE GCMS data
    gc("MHD", "GCMS", "AGAGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)

    # AGAGE Medusa
    gc("MHD", "medusa", "AGAGE", input_directory = input_directory, output_directory = output_directory, version = version, date_range = date_range)


    
def array_job(array_index):
    '''
    Run processing scripts through an array job:
    
    qsub -J 1-50 -k oe -j oe -Wsandbox=PRIVATE process_gcwerks_array.sh
    
    Note that will need to increase the max index to more than 40 if we add more stations/instruments
    
    Check job status using qstat -t <JOBID>[] (note the square brackets)
    '''
    
    
    def wrapper(func, args):
        func(*args)
    
    instrument = [ #AGAGE Medusa
        [gc, ("MHD", "medusa", "AGAGE")],
        [gc, ("CGO", "medusa", "AGAGE")],
        [gc, ("GSN", "medusa", "AGAGE")],
        [gc, ("SDZ", "medusa", "AGAGE")],
        [gc, ("THD", "medusa", "AGAGE")],
        [gc, ("RPB", "medusa", "AGAGE")],
        [gc, ("SMO", "medusa", "AGAGE")],
        [gc, ("SIO", "medusa", "AGAGE")],
        [gc, ("JFJ", "medusa", "AGAGE")],
        [gc, ("CMN", "medusa", "AGAGE")],
        [gc, ("ZEP", "medusa", "AGAGE")],
        # AGAGE GC data
        [gc, ("RPB", "GCMD", "AGAGE")],
        [gc, ("CGO", "GCMD", "AGAGE")],
        [gc, ("MHD", "GCMD", "AGAGE")],
        [gc, ("SMO", "GCMD", "AGAGE")],
        [gc, ("THD", "GCMD", "AGAGE")],
        # AGAGE GCMS data
        [gc, ("CGO", "GCMS", "AGAGE")],
        [gc, ("MHD", "GCMS", "AGAGE")],
        [gc, ("RPB", "GCMS", "AGAGE")],
        [gc, ("SMO", "GCMS", "AGAGE")],
        [gc, ("THD", "GCMS", "AGAGE")],
        [gc, ("JFJ", "GCMS", "AGAGE")],
        [gc, ("CMN", "GCMS", "AGAGE")],
        [gc, ("ZEP", "GCMS", "AGAGE")],
        # AGAGE CRDS data
        [crds, ("RPB", "AGAGE")],
        # GAUGE CRDS data
        [crds, ("HFD", "DECC")],
        [crds, ("BSD", "DECC")],
        # GAUGE GC data
        [gc, ("BSD", "GCMD", "DECC")],
        [gc, ("HFD", "GCMD", "DECC")],
        # DECC CRDS data
        [crds, ("TTA", "DECC")],
        [crds, ("RGL", "DECC")],
        [crds, ("TAC", "DECC")],
        # DECC GC data
        [gc, ("BSD", "GCMD", "DECC")],
        [gc, ("TAC", "GCMD", "DECC")],
        [gc, ("RGL", "GCMD", "DECC")],
        # DECC Medusa
        [gc, ("TAC", "medusa", "DECC")],
        # Bristol CRDS
#        [crds, ("BRI", "DECC")],
        # ICOS
#        [icos, ("TTA", "DECC")],
        [icos, ("MHD", "ICOS")]]
    
    # Return if index is too large for the above list
    if array_index > len(instrument):
        return 0
    
    # Run the relevant script for each station and instrument
    wrapper(instrument[array_index-1][0], instrument[array_index-1][1])

    # Cleanup for particular station (cleans up everything, not just most recently processed instrument)
    cleanup(instrument[array_index-1][1][0])

    
if __name__ == "__main__":

    # AGAGE Medusa
    gc("MHD", "medusa", "AGAGE")
    gc("CGO", "medusa", "AGAGE")
    gc("GSN", "medusa", "AGAGE")
    gc("SDZ", "medusa", "AGAGE")
    gc("THD", "medusa", "AGAGE")
    gc("RPB", "medusa", "AGAGE")
    gc("SMO", "medusa", "AGAGE")
    gc("SIO", "medusa", "AGAGE")
    gc("JFJ", "medusa", "AGAGE")
    gc("CMN", "medusa", "AGAGE")
    gc("ZEP", "medusa", "AGAGE")

    # AGAGE GC data
    gc("RPB", "GCMD", "AGAGE")
    gc("CGO", "GCMD", "AGAGE")
    gc("MHD", "GCMD", "AGAGE")
    gc("SMO", "GCMD", "AGAGE")
    gc("THD", "GCMD", "AGAGE")

    # AGAGE GCMS data
    gc("CGO", "GCMS", "AGAGE")
    gc("MHD", "GCMS", "AGAGE")
    gc("RPB", "GCMS", "AGAGE")
    gc("SMO", "GCMS", "AGAGE")
    gc("THD", "GCMS", "AGAGE")
    gc("JFJ", "GCMS", "AGAGE")
    gc("CMN", "GCMS", "AGAGE")
    gc("ZEP", "GCMS", "AGAGE")

    # AGAGE CRDS data
    crds("RPB", "AGAGE")

    # GAUGE CRDS data
    crds("HFD", "DECC")
    crds("BSD", "DECC")

    # GAUGE GC data
    gc("BSD", "GCMD", "DECC")
    gc("HFD", "GCMD", "DECC")

    # DECC CRDS data
    crds("TTA", "DECC")
    crds("RGL", "DECC")
    crds("TAC", "DECC")

    # DECC GC data
    gc("TAC", "GCMD", "DECC")
    gc("RGL", "GCMD", "DECC")

    # DECC Medusa
    gc("TAC", "medusa", "DECC")

    # ICOS
    icos("TTA", network = "DECC")
    icos("MHD", network = "ICOS")

    cleanup("CGO")
    cleanup("MHD")
    cleanup("RPB")
    cleanup("THD")
    cleanup("SMO")
    cleanup("GSN")
    cleanup("SDZ")
    cleanup("JFJ")
    cleanup("CMN")
    cleanup("ZEP")

    cleanup("TAC")
    cleanup("RGL")
    cleanup("HFD")
    cleanup("BSD")
    cleanup("TTA")
