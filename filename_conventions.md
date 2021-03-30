# Filename conventions for the ACRG repo

## Observation files

The default central location for storing observation files is <data_path>/obs. You can also create your own data directory, if you want to work with non-standard observation data (e.g. pseudo-data).

Within this folder, obs files should be stored in a folder with the site code. E.g. Mace Head data will be in: <data_path>/obs/MHD/

Files must then be named as follows:

```NETWORK-INSTRUMENT_SITECODE_STARTDATE_SPECIES-{INLET-}VERSION.nc```

**NOTE: DO NOT include underscores or hyphens in any of the individual entries**:
- Network: The measurement network, or institution of the data owner (e.g. AGAGE, UoB, etc.)
- Instrument: The measurement instrument. Make sure this is descriptive enough to allow for an instrument change that measures the same gas at the same site (e.g. GCECD, picarro5310)
- Sitecode: 

 - For measurement sites, this is the three-letter site code (typically from the registered code at the GAWSIS data centre). E.g. MHD, CGO, etc. 
 - For satellite data, this is an indication of both the satellite indicator (e.g. gosat, tropomi) and a name for the region the observation files are related to (seperated by a '-') e.g. gosat-india, tropomi-brazil
 
- Startdate: YYYYMMDD of the first data point
- Species: Gas species (e.g. ch4, co2, cfc11). *Don't* include hyphens in here (e.g. cfc11, not cfc-11). *Do* check that species are named consistently with the acrg_species_info.json keys.
- Inlet (**optional**): If multiple inlets are available at a site, use this entry to specify the inlet name, typically the height with an "m". E.g. "100m".
- Version: A version name/number, typically the date that the file was created.

It's important to note that the final group of parameters (```SPECIES-{INLET-}VERSION```), can be either two elements long (no inlet), or three elements (inlet specified).

Examples:

- Mace Head SF6. Mace Head only has one inlet:

```AGAGE-GCMSMedusa_MHD_20031115_sf6-20201204.nc```

- Bilsdale CH4 on the Picarro, at the 248m inlet:

```DECC-picarro_BSD_20140130_ch4-248m-20201204.nc```

## Footprint files

Default path: <data_path>/LPDM/fp_NAME/ or <data_path>/LPDM/fp_NAME_high_time_res/

Footprint file names are defined in acrg_name/name.py::filenames(...) as

```[fp_directory]/domain/site*-height-species*domain*ym*.nc ``` or ```[fp_directory]/domain/site*-height_domain*ym*.nc ```

- site: sitecode
- height: release height of footprint (may not be same as actual inlet height, check site json)
- species: (optional) for short-lived species to account for decay during 30 day footprint, leave out for all long lived tracer species
- domain: NAME domain for file
- ym: yearmonth format date of the footprint file

## Boundary condition files

Default path: <data_path>/LPDM/bc/

```[bc_directory]/domain/species_*.nc ```

- domain: NAME domain
- species: (all lowercase) species name

## Basis files

Default path: <data_path>/LPDM/basis_function/ or <data_path>/LPDM/bc_basis_function/ for bc basis

```[basis_directory]/domain/basis_domain*.nc ```

- domain: NAME domain
- basis: the basis case to be used, eg quadtree

## Emissions files

Default path: <data_path>/LPDM/emissions/

```[emissions_directory]/domain/species_*.nc ```

- domain: NAME domain
- species: (all lowercase) species to grab emissions for. Also used to signify sectors etc
