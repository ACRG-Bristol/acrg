# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]
Use this section to keep track of any intermediate changes that can then be moved to a versioned section.

Updated pinned packages, now expects Python 3.9

### Added

### Changed

### Removed


## [0.3.0] - 2022-04-01
### Added
- hbmcmc code now has a function (and relevant additions to output) that allows the inversion to be rerun (i.e. reproduced)
  using only the output as inputs and ACRG repository.
- Country mask code has now been updated to allow ocean territories (Economic Exclusion Zone, EEZ) to be included.
- hbmcmc code now has options to use HiTRes processes for the set up
- Add function for embedding regional fields into larger fields, in acrg.name.emissions_helper_func

### Changed
- directory structure of acrg package has changed. Many import statements will likely need modifying
- the filename structure of the footprint nc files created in acrg_name.process has changed and the new filenames are detailed in filename_convention.md
- the footprint files now have an explicit convention for labeling a met_model such as 'UKV'
- allows species lifetime to be specified monthly or annually in species_info.json
- allows CO2 footprint files in fp_NAME_high_time_res to be processed using the acrg_name.process code by calling species = 'CO2'
- updated acrg_hbmcmc files to accommodate new footprint naming structure and to read in met_model through the .ini file
- no changes will need to be made to user workflow
- Some additions to hbmcmc output and changed string format for some attributes
- update name.name.footprints_data_merge to accommodate high time res emissions, e.g. CO2.
- NAME process script has been rewritten to not rely on load_NAME proprietry code

### Removed
- removed GCWerks modules


## [0.2.0] - 2020-04-01
### Added
- Ability to convert calibration scale in get_obs
- New "defaults" file that specifies inlets and instruments to use for particular time periods
- An obs.db SQLite database that specifies the location of all obs files and basic details about their contents (species, inlet, time range, etc.)
- notebooks directory for Jupyter notebooks
- notebooks/tutorials directory for notebook based tutorials
- a tmp directory to store random job script output files
- added a dev environment that includes spyder and a lighter environment that does not
- acrg_BC for creating boundary conditions
- 2021/06/01, LMW: Added ability to have an offset between sites in hbmcmc code. 


### Changed
- get_single_site now returns a list of xarray datasets, one for each combination of inlet and site. If defaults are specified, the list will contain the default instruments and inlets for each period
- get_obs now returns a dictionary containing lists of datasets
- calibration scale and inlet are now attributes to obs datasets (e.g. ds.attrs["scale"])
- fp_data_merge now works with new get_obs object
- The flux function will now look for species-total_*.nc named files first and then look for species_*.nc files. This will not be able to read both files. This can still accept an more explicit source such as co2-ff_*.nc as an alternative to this. 
- arviz package version pinned to prevent conflict with pymc3 version

### Removed
- N/A
