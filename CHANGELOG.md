# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]
Use this section to keep track of any intermediate changes that can then be moved to a versioned section.

## [0.2.0] - 2020-XX-XX
### Added
- Ability to convert calibration scale in get_obs
- New "defaults" file that specifies inlets and instruments to use for particular time periods
- An obs.db SQLite database that specifies the location of all obs files and basic details about their contents (species, inlet, time range, etc.)
- notebooks directory for Jupyter notebooks
- notebooks/tutorials directory for notebook based tutorials
- a tmp directory to store random job script output files

### Changed
- get_single_site now returns a list of xarray datasets, one for each combination of inlet and site. If defaults are specified, the list will contain the default instruments and inlets for each period
- get_obs now returns a dictionary containing lists of datasets
- calibration scale and inlet are now attributes to obs datasets (e.g. ds.attrs["scale"])
- fp_data_merge now works with new get_obs object

### Removed
- N/A