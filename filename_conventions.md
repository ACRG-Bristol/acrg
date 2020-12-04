# Filename conventions for the ACRG repo

## Observation files

The central location for storing observation files on BP1 is /work/chxmr/shared/obs. You can also create your own data directory, if you want to work with non-standard observation data (e.g. pseudo-data).

Within this folder, obs files should be stored in a folder with the site code. E.g. Mace Head data will be in:

/work/chxmr/shared/obs/MHD/

Files must then be named as follows:

```NETWORK-INSTRUMENT_SITECODE_STARTDATE_SPECIES-{INLET-}VERSION.nc```

- Network: The measurement network, or institution of the data owner (e.g. AGAGE, UoB, etc.)
- Instrument: The measurement instrument. Make sure this is descriptive enough to allow for an instrument change that measures the same gas at the same site (e.g. GCECD, picarro5310)
- Sitecode: Three-letter site code (typically from the registered code at the GAWSIS data centre). E.g. MHD, CGO, etc.
- Startdate: YYYYMMDD of the first data point
- Species: Gas species (e.g. ch4, co2, cfc11). *Don't* include hyphens in here (e.g. cfc11, not cfc-11)
- Inlet (**optional**): If multiple inlets are available at a site, use this entry to specify the inlet name, typically the height with an "m". E.g. "100m".
- Version: A version name/number, typically the date that the file was created.

It's important to note that the final group of parameters, separated by "-", can be either two elements long (no inlet), or three elements (inlet specified).

Examples:

- Mace Head SF6. Mace Head only has one inlet:

```AGAGE-GCMSMedusa_MHD_20031115_sf6-20201204.nc```

- Bilsdale CH4 on the Picarro, at the 248m inlet:

```DECC-picarro_BSD_20140130_ch4-248m-20201204.nc```

## Footprint files



## Boundary condition files



## Emissions files


